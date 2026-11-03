"""Environments (media, oxygen, nitrogen source, ...) for yeast-GEM.

Each environment is a file ``data/conditions/<name>.yml``. The MATLAB
``applyEnvironment`` reads the same files and applies the same steps in the
same order, so both languages give identical models; the numbers in the files
are read the same way in both (``1e3`` is 1000, ``010`` is 10, ``yes`` is a
string).

This module needs only cobrapy and pyyaml, and from ``yeastgem`` only
``yeastgem.paths`` (standard library), so applying environments does not
depend on raven-toolbox: importing ``yeastgem.conditions`` does not import it.
The environment is checked before the model is changed: an invalid
environment, including an identifier that is not in the model, raises
ValueError and leaves the model as it was.

Steps, in this order:

``amino_acid_ratio``
    ``aerobic`` or ``anaerobic`` column of
    ``data/physiology/aminoAcid_Bjorkeroth2020.tsv``; the protein
    pseudoreaction is rescaled to its previous mass (``changeAminoAcidRatio.m``).
``prelude.reset_exchanges``
    ``out``, ``in`` or ``all``: set those exchange reactions to (0, 1000), which
    closes their uptake and opens their secretion. ``out`` are the reactions
    without products (``A =>``), ``in`` those without substrates (``=> A``) and
    ``all`` both, as RAVEN's getExchangeRxns.
``cofactor_pseudoreaction``
    remove metabolites and recompute the charge-balancing H+ (a metabolite
    without charge counts as 0, with a warning).
``biomass_stoichiometry_delta``
    add coefficients to a reaction.
``bounds``
    set lb and/or ub of listed reactions. Entries for one reaction are
    combined in order, the later value of a key winning.
``expected_uptake_count``
    warn if a different number of lb -1000 bounds was applied.
"""
from __future__ import annotations

import csv
import functools
import re
import warnings
from pathlib import Path
from typing import Any

import cobra
import yaml

from yeastgem.paths import find_repo_root

_PROTEIN_MET_NAME = "protein"  # product of the protein pseudoreaction

# Atomic weights as in sumBioMass.m; the generic residue R weighs 0.
_ELEMENT_WEIGHTS = {
    "C": 12.01, "H": 1.008, "N": 14.007, "O": 15.999, "P": 30.974,
    "S": 32.06, "R": 0.0, "Fe": 55.845, "K": 39.098, "Na": 22.99,
    "Cl": 35.45, "Mn": 54.938, "Zn": 65.38, "Ca": 40.078, "Mg": 24.305,
    "Cu": 63.546,
}
_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]*)(\d*)")

# Scalars as the MATLAB reader reads them: every number is a decimal number
# (1e3, 010), only true and false are booleans, yes and no stay strings.
_NUMBER = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")
_BOOLEAN = re.compile(r"^(?:true|True|TRUE|false|False|FALSE)$")
_REPLACED_TAGS = {f"tag:yaml.org,2002:{t}" for t in ("bool", "float", "int", "timestamp")}


class _Loader(yaml.SafeLoader):
    """Reads data/conditions the way the MATLAB reader does."""


_Loader.yaml_implicit_resolvers = {
    first: [(tag, regexp) for tag, regexp in resolvers if tag not in _REPLACED_TAGS]
    for first, resolvers in yaml.SafeLoader.yaml_implicit_resolvers.items()
}
_Loader.add_implicit_resolver("tag:yaml.org,2002:float", _NUMBER, list("-+0123456789."))
_Loader.add_implicit_resolver("tag:yaml.org,2002:bool", _BOOLEAN, list("tTfF"))


@functools.lru_cache(maxsize=1)
def _ids() -> dict[str, Any]:
    """The canonical identifiers of ``data/yeastgem/ids.yml``."""
    with open(find_repo_root() / "data" / "yeastgem" / "ids.yml", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def load_condition(name: str, *, conditions_dir: Path | None = None) -> dict[str, Any]:
    """Read ``data/conditions/<name>.yml`` (or a path to such a file)."""
    path = Path(name)
    if path.suffix != ".yml":
        base = Path(conditions_dir) if conditions_dir else find_repo_root() / "data" / "conditions"
        path = base / f"{name}.yml"
    if not path.is_file():
        raise FileNotFoundError(f"No environment file {path}")
    with open(path, encoding="utf-8") as handle:
        return yaml.load(handle, Loader=_Loader) or {}


def apply(model: cobra.Model, name: str) -> cobra.Model:
    """Apply the environment ``name`` to ``model`` in place and return it."""
    return apply_condition(model, load_condition(name))


def apply_condition(model: cobra.Model, cfg: dict[str, Any]) -> cobra.Model:
    """Apply a parsed environment (a dict as read by :func:`load_condition`).

    Every identifier is looked up and every value checked before the first
    change, so an invalid environment raises ValueError with the model as it was.
    """
    reset = (cfg.get("prelude") or {}).get("reset_exchanges")
    if reset is not None:
        if not isinstance(reset, str) or reset.lower() not in ("in", "out", "all"):
            raise ValueError(f"reset_exchanges must be in, out or all, not {reset!r}")
        reset = reset.lower()
    ratio = cfg.get("amino_acid_ratio")
    if ratio is not None and ratio not in ("aerobic", "anaerobic"):
        raise ValueError(f"amino_acid_ratio must be aerobic or anaerobic, not {ratio!r}")

    exchange = {}
    if reset:
        for rxn in model.reactions:
            no_products = not any(c > 0 for c in rxn.metabolites.values())
            no_substrates = not any(c < 0 for c in rxn.metabolites.values())
            if {"out": no_products, "in": no_substrates}.get(reset, no_products or no_substrates):
                exchange[rxn.id] = rxn
    bounds = _checked_bounds(model, cfg.get("bounds") or [], exchange)
    amino_acids = _amino_acid_plan(model, aerobic=ratio == "aerobic") if ratio else None

    cp = cfg.get("cofactor_pseudoreaction")
    if cp:
        cp_rxn = _find(model.reactions, cp.get("rxn_id"), "Reaction")
        cp_remove = [_find(model.metabolites, e.get("met"), "Metabolite")
                     for e in cp.get("remove_mets") or []]
        cp_balance = (_find(model.metabolites, cp["charge_balance_met"], "Metabolite")
                      if cp.get("charge_balance_met") else None)
    delta = cfg.get("biomass_stoichiometry_delta")
    if delta:
        delta_rxn = _find(model.reactions, delta.get("rxn_id"), "Reaction")
        delta_add = {_find(model.metabolites, e.get("met"), "Metabolite"): float(e["coef"])
                     for e in delta.get("add") or []}

    # Nothing above changed the model; from here on it is changed.
    if amino_acids:
        _apply_amino_acid_plan(model, *amino_acids)

    for rxn in exchange.values():
        rxn.bounds = (0.0, 1000.0)

    if cp:
        for met in cp_remove:
            _set_coefficient(cp_rxn, met, 0.0)
        if cp_balance is not None:
            _rebalance_charge(cp_rxn, cp_balance)

    if delta:
        delta_rxn.add_metabolites(delta_add, combine=True)

    n_uptake = 0
    for rxn, lb, ub, has_lb in bounds:
        rxn.bounds = (lb, ub)
        n_uptake += has_lb and lb == -1000
    expected = cfg.get("expected_uptake_count")
    if expected is not None and n_uptake != expected:
        warnings.warn(f"Expected {expected} uptake reactions, applied {n_uptake}.", stacklevel=2)
    return model


def _find(collection, identifier: Any, kind: str):
    """The reaction or metabolite ``identifier``; ValueError if it is not in the model."""
    if not isinstance(identifier, str) or identifier not in collection:
        raise ValueError(f"{kind} {identifier!r} is not in the model.")
    return collection.get_by_id(identifier)


def _checked_bounds(
    model: cobra.Model, entries: list[dict[str, Any]], exchange: dict[str, cobra.Reaction]
) -> list[tuple[cobra.Reaction, float, float, bool]]:
    """The bounds to set, one per reaction, as (reaction, lb, ub, lb given).

    Entries for one reaction are combined in order, starting from the bounds
    after the exchange reset. Raises ValueError for a reaction that is not in
    the model, a non-numeric bound or lb > ub."""
    combined: dict[str, list] = {}
    for entry in entries:
        rxn = _find(model.reactions, entry.get("rxn"), "Reaction")
        if rxn.id not in combined:
            lb, ub = (0.0, 1000.0) if rxn.id in exchange else rxn.bounds
            combined[rxn.id] = [rxn, lb, ub, False]
        item = combined[rxn.id]
        for position, key in ((1, "lb"), (2, "ub")):
            if key in entry:
                value = entry[key]
                if isinstance(value, bool) or not isinstance(value, (int, float)) \
                        or value != value or abs(value) == float("inf"):
                    raise ValueError(f"{rxn.id}: {key} must be a finite number, not {value!r}")
                item[position] = float(value)
                item[3] = item[3] or key == "lb"
    for rxn, lb, ub, _ in combined.values():
        if lb > ub:
            raise ValueError(f"{rxn.id}: lb {lb} is larger than ub {ub}")
    return [tuple(item) for item in combined.values()]


def change_amino_acid_ratio(
    model: cobra.Model, *, aerobic: bool = True, aa_tsv: Path | str | None = None
) -> cobra.Model:
    """Replace the amino-acid ratios of the protein pseudoreaction and rescale
    it to the previous protein mass (``changeAminoAcidRatio.m``)."""
    _apply_amino_acid_plan(model, *_amino_acid_plan(model, aerobic=aerobic, aa_tsv=aa_tsv))
    return model


def _amino_acid_plan(
    model: cobra.Model, *, aerobic: bool, aa_tsv: Path | str | None = None
) -> tuple[cobra.Reaction, list[tuple[cobra.Metabolite, cobra.Metabolite, float]],
           float, cobra.Metabolite]:
    """What :func:`change_amino_acid_ratio` sets, with every metabolite looked up
    and weighed, so that a missing metabolite or formula raises before the model
    is changed."""
    ids = _ids()
    rxn = _find(model.reactions, ids["protein_rxn"], "Reaction")
    proton = _find(model.metabolites, ids["proton_met"], "Metabolite")
    column = 4 if aerobic else 5
    path = Path(aa_tsv) if aa_tsv else (
        find_repo_root() / "data" / "physiology" / "aminoAcid_Bjorkeroth2020.tsv")
    ratios = []
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for row in reader:
            if row and row[0].strip():
                ratios.append((_find(model.metabolites, row[1], "Metabolite"),
                               _find(model.metabolites, row[2], "Metabolite"),
                               float(row[column])))
    for substrate, _, _ in ratios:
        _formula_weight(substrate)
    return rxn, ratios, _protein_mass(model), proton


def _apply_amino_acid_plan(model, rxn, ratios, target, proton) -> None:
    for substrate, product, ratio in ratios:
        _set_coefficient(rxn, substrate, -ratio)
        _set_coefficient(rxn, product, ratio)
    factor = target / _protein_mass(model)
    deltas = {
        m: (factor - 1.0) * c for m, c in rxn.metabolites.items() if m.name != _PROTEIN_MET_NAME
    }
    if deltas:
        rxn.add_metabolites(deltas, combine=True)
    _rebalance_charge(rxn, proton, unknown_as_zero=True)


# --- helpers ----------------------------------------------------------

def _set_coefficient(rxn: cobra.Reaction, met: cobra.Metabolite, value: float) -> None:
    """Set an absolute coefficient (MATLAB ``S(i,j) = value``); 0 removes ``met``."""
    change = float(value) - rxn.metabolites.get(met, 0.0)
    if change != 0.0:
        rxn.add_metabolites({met: change}, combine=True)


def _rebalance_charge(
    rxn: cobra.Reaction, balance_met: cobra.Metabolite, *, unknown_as_zero: bool = False
) -> None:
    """Set ``balance_met``'s coefficient so that ``rxn`` is charge balanced.
    A metabolite without charge counts as 0 ('omitnan' in MATLAB), with a
    warning unless ``unknown_as_zero``."""
    _set_coefficient(rxn, balance_met, 0.0)
    unknown = sorted(m.id for m in rxn.metabolites if m.charge is None)
    if unknown and not unknown_as_zero:
        warnings.warn(
            f"Charge balance of {rxn.id}: no charge for {', '.join(unknown)}, counted as 0.",
            stacklevel=3,
        )
    total = sum((m.charge or 0) * c for m, c in rxn.metabolites.items())
    _set_coefficient(rxn, balance_met, -total)


def _protein_mass(model: cobra.Model) -> float:
    """Protein fraction [g/gDW] (sumBioMass.m): substrate MWs of the protein
    pseudoreaction, two protons (2.016) removed per charged tRNA."""
    rxn = model.reactions.get_by_id(_ids()["protein_rxn"])
    return sum(
        -c * (_formula_weight(m) - 2.016) for m, c in rxn.metabolites.items() if c < 0
    ) / 1000.0


def _formula_weight(met: cobra.Metabolite) -> float:
    if not met.formula:
        raise ValueError(f"Biomass metabolite {met.id} has no formula.")
    weight = 0.0
    for element, count in _FORMULA_TOKEN.findall(met.formula):
        if element not in _ELEMENT_WEIGHTS:
            raise ValueError(f"Unknown element {element!r} in formula {met.formula!r}.")
        weight += (int(count) if count else 1) * _ELEMENT_WEIGHTS[element]
    return weight
