"""Apply named condition presets to a cobra.Model.

Mirrors `code/applyCondition.m`; the YAML files in `data/conditions/`
are the single source of truth shared with the MATLAB loader.

Phase-2 scope: applies `prelude` (reset exchanges), `cofactor_pseudoreaction`
(remove mets + H+ charge balance), `biomass_stoichiometry_delta`, and the
per-reaction `bounds` list. The `amino_acid_ratio` step depends on the
biomass machinery and is deferred to Tier 2 — calling
``apply(model, 'anaerobic')`` therefore raises `NotImplementedError`
until phase 4 lands. Apply the anaerobic condition via MATLAB for now.
"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import cobra
import yaml

from yeastgem.io import REPO_PATH

_CONDITIONS_DIR = REPO_PATH / "data" / "conditions"


def load_condition(name: str, *, conditions_dir: Path | None = None) -> dict[str, Any]:
    """Load `data/conditions/<name>.yml` as a plain dict.

    Used by both `apply` and downstream callers that want to inspect the
    raw condition spec.
    """
    base = conditions_dir or _CONDITIONS_DIR
    path = base / f"{name}.yml"
    if not path.exists():
        raise FileNotFoundError(f"No such condition: {name} (looked for {path})")
    with open(path) as f:
        return yaml.safe_load(f)


def apply(model: cobra.Model, name: str) -> cobra.Model:
    """Apply the named condition to ``model`` in-place and return it.

    Mirrors ``applyCondition.m``. Raises ``NotImplementedError`` if the
    condition uses ``amino_acid_ratio`` (Tier 2 dependency).
    """
    cfg = load_condition(name)

    _apply_prelude(model, cfg)
    _apply_cofactor_pseudoreaction(model, cfg)
    _check_amino_acid_ratio_unsupported(cfg, name)
    _apply_biomass_stoichiometry_delta(model, cfg)
    n_uptake = _apply_bounds(model, cfg)
    _check_uptake_count(cfg, n_uptake)

    return model


# --- step implementations ---------------------------------------------

def _apply_prelude(model: cobra.Model, cfg: dict[str, Any]) -> None:
    prelude = cfg.get("prelude")
    if not prelude:
        return
    if prelude.get("reset_exchanges"):
        # cobrapy does not distinguish RAVEN's 'in'/'out' direction; all
        # exchanges are reset to (lb=0, ub=1000). For yeast-GEM the two
        # produce the same bound set.
        for rxn in model.exchanges:
            rxn.lower_bound = 0
            rxn.upper_bound = 1000


def _apply_cofactor_pseudoreaction(model: cobra.Model, cfg: dict[str, Any]) -> None:
    cp = cfg.get("cofactor_pseudoreaction")
    if not cp:
        return
    rxn = model.reactions.get_by_id(cp["rxn_id"])
    for entry in cp.get("remove_mets", []):
        met = model.metabolites.get_by_id(entry["met"])
        _set_coefficient(rxn, met, 0.0)
    balance_met_id = cp.get("charge_balance_met")
    if balance_met_id:
        balance_met = model.metabolites.get_by_id(balance_met_id)
        _set_coefficient(rxn, balance_met, 0.0)
        total_charge = sum(
            (m.charge or 0) * coef
            for m, coef in rxn.metabolites.items()
        )
        _set_coefficient(rxn, balance_met, -total_charge)


def _check_amino_acid_ratio_unsupported(cfg: dict[str, Any], name: str) -> None:
    if "amino_acid_ratio" in cfg:
        raise NotImplementedError(
            f"Condition {name!r} uses `amino_acid_ratio`; the Python "
            "implementation depends on the biomass module and is deferred "
            "to Tier 2 (see code/python/PORTING_PLAN.md). Apply this "
            "condition via MATLAB until then."
        )


def _apply_biomass_stoichiometry_delta(model: cobra.Model, cfg: dict[str, Any]) -> None:
    delta = cfg.get("biomass_stoichiometry_delta")
    if not delta:
        return
    rxn = model.reactions.get_by_id(delta["rxn_id"])
    for entry in delta.get("add", []):
        met = model.metabolites.get_by_id(entry["met"])
        rxn.add_metabolites({met: float(entry["coef"])}, combine=True)


def _apply_bounds(model: cobra.Model, cfg: dict[str, Any]) -> int:
    n_uptake = 0
    for entry in cfg.get("bounds", []):
        try:
            rxn = model.reactions.get_by_id(entry["rxn"])
        except KeyError:
            warnings.warn(
                f"Reaction {entry['rxn']!r} not found in model; skipping.",
                stacklevel=3,
            )
            continue
        new_lb = float(entry["lb"]) if "lb" in entry else rxn.lower_bound
        new_ub = float(entry["ub"]) if "ub" in entry else rxn.upper_bound
        _set_bounds(rxn, new_lb, new_ub)
        # `expected_uptake_count` mirrors the MATLAB minimal_Y6 sanity
        # check, which counted the 15 desiredExchanges (all set to
        # lb=-1000). Glucose (lb=-1) is not part of that group.
        if entry.get("lb") == -1000:
            n_uptake += 1
    return n_uptake


def _set_bounds(rxn: cobra.Reaction, lb: float, ub: float) -> None:
    """Set both bounds, bypassing cobra's lb<=ub validator if necessary.

    Some legacy conditions (glycine_nitrogen, nitrogen_limitation) set
    lb=1000, ub=0 on the glycine cleavage reactions. The resulting model
    is infeasible, but reproducing it exactly is required for the
    lock-step parity check between MATLAB and Python.

    cobra's `lower_bound`/`upper_bound` setters and the `bounds`
    setter both reject ``lb > ub``. The only way to land on a matching
    state is to write the underlying private attributes directly.
    """
    if lb > ub:
        rxn._lower_bound = lb
        rxn._upper_bound = ub
    else:
        rxn.bounds = (lb, ub)


def _check_uptake_count(cfg: dict[str, Any], n_uptake: int) -> None:
    expected = cfg.get("expected_uptake_count")
    if expected is None:
        return
    if n_uptake != expected:
        warnings.warn(
            f"Expected {expected} uptake reactions, applied {n_uptake}. "
            "Some referenced reactions may be missing from the model.",
            stacklevel=3,
        )


# --- helpers ----------------------------------------------------------

def _set_coefficient(rxn: cobra.Reaction, met: cobra.Metabolite, value: float) -> None:
    """Set the stoichiometric coefficient of ``met`` in ``rxn`` to ``value``.

    Mirrors ``model.S(metIdx, rxnIdx) = value`` in MATLAB. cobra's
    ``add_metabolites`` is the only public mutation point; we go via
    `combine=True` with `(value - current)` to land on the desired value
    without depending on cobra's removal-on-zero behaviour.
    """
    current = rxn.metabolites.get(met, 0.0)
    delta = float(value) - current
    if delta != 0:
        rxn.add_metabolites({met: delta}, combine=True)
