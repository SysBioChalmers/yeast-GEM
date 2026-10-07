"""Environments (media, oxygen, nitrogen source, ...) for yeast-GEM.

Each environment is a file ``data/conditions/<name>.yml``. The MATLAB
``applyEnvironment`` reads the same files and applies the same steps in the
same order, so both languages give identical models.

This module needs only cobrapy and pyyaml, and from ``yeastgem`` only
``yeastgem.paths`` (standard library), so applying environments does not
depend on raven-toolbox. The environment is checked before the model is
changed: an invalid environment raises ValueError and leaves the model as it
was.

Steps, in this order:

``amino_acid_ratio``
    ``aerobic`` or ``anaerobic`` column of
    ``data/physiology/aminoAcid_Bjorkeroth2020.tsv``; the protein
    pseudoreaction is rescaled to its previous mass (``changeAminoAcidRatio.m``).
``prelude.reset_exchanges``
    ``out``, ``in`` or ``all``: set those exchange reactions (no products,
    no substrates, or either; as RAVEN's getExchangeRxns) to (0, 1000).
``cofactor_pseudoreaction``
    remove metabolites and recompute the charge-balancing H+ (a metabolite
    without charge counts as 0, with a warning).
``biomass_stoichiometry_delta``
    add coefficients to a reaction.
``bounds``
    set lb and/or ub of listed reactions.
``expected_uptake_count``
    warn if a different number of lb -1000 bounds was applied.
"""
from __future__ import annotations

import csv
import re
import warnings
from pathlib import Path
from typing import Any

import cobra
import yaml

from yeastgem.paths import find_repo_root

_PROTON_ID = "s_0794"          # H+ [cytoplasm]
_PROTEIN_RXN = "r_4047"        # protein pseudoreaction
_PROTEIN_MET_NAME = "protein"  # its product

# Atomic weights as in sumBioMass.m; the generic residue R weighs 0.
_ELEMENT_WEIGHTS = {
    "C": 12.01, "H": 1.008, "N": 14.007, "O": 15.999, "P": 30.974,
    "S": 32.06, "R": 0.0, "Fe": 55.845, "K": 39.098, "Na": 22.99,
    "Cl": 35.45, "Mn": 54.938, "Zn": 65.38, "Ca": 40.078, "Mg": 24.305,
    "Cu": 63.546,
}
_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]*)(\d*)")


def load_condition(name: str, *, conditions_dir: Path | None = None) -> dict[str, Any]:
    """Read ``data/conditions/<name>.yml`` (or a path to such a file)."""
    path = Path(name)
    if path.suffix != ".yml":
        base = Path(conditions_dir) if conditions_dir else find_repo_root() / "data" / "conditions"
        path = base / f"{name}.yml"
    if not path.is_file():
        raise FileNotFoundError(f"No environment file {path}")
    with open(path, encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def apply(model: cobra.Model, name: str) -> cobra.Model:
    """Apply the environment ``name`` to ``model`` in place and return it."""
    return apply_condition(model, load_condition(name))


def apply_condition(model: cobra.Model, cfg: dict[str, Any]) -> cobra.Model:
    """Apply a parsed environment (a dict as read by :func:`load_condition`)."""
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

    if ratio is not None:
        change_amino_acid_ratio(model, aerobic=ratio == "aerobic")

    for rxn in exchange.values():
        rxn.bounds = (0.0, 1000.0)

    cp = cfg.get("cofactor_pseudoreaction")
    if cp:
        rxn = model.reactions.get_by_id(cp["rxn_id"])
        for entry in cp.get("remove_mets") or []:
            _set_coefficient(rxn, model.metabolites.get_by_id(entry["met"]), 0.0)
        if cp.get("charge_balance_met"):
            _rebalance_charge(rxn, model.metabolites.get_by_id(cp["charge_balance_met"]))

    delta = cfg.get("biomass_stoichiometry_delta")
    if delta:
        rxn = model.reactions.get_by_id(delta["rxn_id"])
        rxn.add_metabolites(
            {model.metabolites.get_by_id(e["met"]): float(e["coef"])
             for e in delta.get("add") or []},
            combine=True,
        )

    n_uptake = 0
    for rxn, lb, ub, has_lb in bounds:
        rxn.bounds = (lb, ub)
        n_uptake += has_lb and lb == -1000
    expected = cfg.get("expected_uptake_count")
    if expected is not None and n_uptake != expected:
        warnings.warn(f"Expected {expected} uptake reactions, applied {n_uptake}.", stacklevel=2)
    return model


def _checked_bounds(
    model: cobra.Model, entries: list[dict[str, Any]], exchange: dict[str, cobra.Reaction]
) -> list[tuple[cobra.Reaction, float, float, bool]]:
    """The bounds to set, as (reaction, lb, ub, lb given). Raises ValueError
    for a non-numeric bound or lb > ub, taking the exchange reset into account."""
    out = []
    for entry in entries:
        if entry["rxn"] not in model.reactions:
            warnings.warn(f"Reaction {entry['rxn']} not in the model; skipped.", stacklevel=3)
            continue
        rxn = model.reactions.get_by_id(entry["rxn"])
        lb, ub = (0.0, 1000.0) if rxn.id in exchange else rxn.bounds
        for key in ("lb", "ub"):
            if key in entry:
                value = entry[key]
                if isinstance(value, bool) or not isinstance(value, (int, float)) \
                        or value != value or abs(value) == float("inf"):
                    raise ValueError(f"{rxn.id}: {key} must be a finite number, not {value!r}")
        lb = float(entry.get("lb", lb))
        ub = float(entry.get("ub", ub))
        if lb > ub:
            raise ValueError(f"{rxn.id}: lb {lb} is larger than ub {ub}")
        out.append((rxn, lb, ub, "lb" in entry))
    return out


def change_amino_acid_ratio(
    model: cobra.Model, *, aerobic: bool = True, aa_tsv: Path | str | None = None
) -> cobra.Model:
    """Replace the amino-acid ratios of the protein pseudoreaction and rescale
    it to the previous protein mass (``changeAminoAcidRatio.m``)."""
    target = _protein_mass(model)
    rxn = model.reactions.get_by_id(_PROTEIN_RXN)
    column = 4 if aerobic else 5
    path = Path(aa_tsv) if aa_tsv else (
        find_repo_root() / "data" / "physiology" / "aminoAcid_Bjorkeroth2020.tsv")
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for row in reader:
            if row and row[0].strip():
                ratio = float(row[column])
                _set_coefficient(rxn, model.metabolites.get_by_id(row[1]), -ratio)
                _set_coefficient(rxn, model.metabolites.get_by_id(row[2]), ratio)
    factor = target / _protein_mass(model)
    deltas = {
        m: (factor - 1.0) * c for m, c in rxn.metabolites.items() if m.name != _PROTEIN_MET_NAME
    }
    if deltas:
        rxn.add_metabolites(deltas, combine=True)
    _rebalance_charge(rxn, model.metabolites.get_by_id(_PROTON_ID), unknown_as_zero=True)
    return model


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
    rxn = model.reactions.get_by_id(_PROTEIN_RXN)
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
