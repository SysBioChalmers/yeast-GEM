"""Tests for ``yeastgem.conditions`` (data-driven condition presets)."""
from __future__ import annotations

import pytest

from yeastgem import compare_models, conditions

# --- data-file shape checks (no model load needed) --------------------

def test_minimal_Y6_loads():
    cfg = conditions.load_condition("minimal_Y6")
    assert cfg["name"] == "minimal_Y6"
    assert cfg["prelude"]["reset_exchanges"] == "out"
    assert cfg["expected_uptake_count"] == 15


def test_anaerobic_loads():
    cfg = conditions.load_condition("anaerobic")
    assert cfg["name"] == "anaerobic"
    assert cfg["amino_acid_ratio"] == "anaerobic"
    assert cfg["cofactor_pseudoreaction"]["rxn_id"] == "r_4598"
    assert cfg["biomass_stoichiometry_delta"]["rxn_id"] == "r_4041"


def test_glycine_nitrogen_loads():
    cfg = conditions.load_condition("glycine_nitrogen")
    assert cfg["name"] == "glycine_nitrogen"
    assert {b["rxn"] for b in cfg["bounds"]} == {"r_0501", "r_0507", "r_0509"}


def test_nitrogen_limitation_loads():
    cfg = conditions.load_condition("nitrogen_limitation")
    assert cfg["name"] == "nitrogen_limitation"


def test_unknown_condition_raises():
    with pytest.raises(FileNotFoundError):
        conditions.load_condition("does_not_exist")


# --- application checks (need the model) ------------------------------

def test_apply_glycine_nitrogen_sets_bounds(model):
    mutated = model.copy()
    conditions.apply(mutated, "glycine_nitrogen")
    for rxn_id in ("r_0501", "r_0507", "r_0509"):
        rxn = mutated.reactions.get_by_id(rxn_id)
        assert rxn.lower_bound == 1000
        assert rxn.upper_bound == 0


def test_apply_nitrogen_limitation_sets_bounds(model):
    mutated = model.copy()
    conditions.apply(mutated, "nitrogen_limitation")
    assert mutated.reactions.get_by_id("r_0472").upper_bound == 1000
    for rxn_id in ("r_0501", "r_0507", "r_0509"):
        assert mutated.reactions.get_by_id(rxn_id).lower_bound == 1000


def test_apply_minimal_Y6_caps_glucose_and_zeros_bicarbonate(model):
    mutated = model.copy()
    conditions.apply(mutated, "minimal_Y6")
    glucose = mutated.reactions.get_by_id("r_1714")
    assert glucose.lower_bound == -1
    bicarbonate = mutated.reactions.get_by_id("r_1663")
    assert bicarbonate.lower_bound == 0
    assert bicarbonate.upper_bound == 0
    # Allowed uptakes (sample a few)
    for rxn_id in ("r_1654", "r_1992", "r_2005", "r_2060"):
        assert mutated.reactions.get_by_id(rxn_id).lower_bound == -1000


def test_apply_minimal_Y6_resets_all_exchanges(model):
    """Prelude should set all "out" exchanges to (lb=0, ub=1000) before
    the targeted overrides — confirmed by the bicarbonate/oxygen path."""
    mutated = model.copy()
    # Pick an exchange that the condition does NOT touch and verify the
    # prelude zeroed its lb (uptake blocked) and capped its ub at 1000.
    untouched_exchange = next(
        r for r in mutated.exchanges
        if r.id not in {b["rxn"] for b in conditions.load_condition("minimal_Y6")["bounds"]}
    )
    conditions.apply(mutated, "minimal_Y6")
    assert untouched_exchange.lower_bound == 0
    assert untouched_exchange.upper_bound == 1000


def test_apply_anaerobic_raises_until_tier2(model):
    """Anaerobic uses amino_acid_ratio; Python implementation deferred to
    Tier 2. The error message must point to the plan."""
    mutated = model.copy()
    with pytest.raises(NotImplementedError, match="amino_acid_ratio"):
        conditions.apply(mutated, "anaerobic")


def test_apply_is_idempotent_for_glycine(model):
    """Applying glycine_nitrogen twice must produce the same model."""
    once = model.copy()
    conditions.apply(once, "glycine_nitrogen")
    twice = model.copy()
    conditions.apply(twice, "glycine_nitrogen")
    conditions.apply(twice, "glycine_nitrogen")
    report = compare_models(once, twice)
    assert report.equal, report


# --- partial-anaerobic check (everything except amino_acid_ratio) -----

def test_anaerobic_cofactor_delta_changes_S_matrix(model):
    """Manually invoke the parts of anaerobic Python supports and verify
    the cofactor pseudoreaction loses heme a (s_3714)."""
    mutated = model.copy()
    cofac = mutated.reactions.get_by_id("r_4598")
    heme = mutated.metabolites.get_by_id("s_3714")
    assert heme in cofac.metabolites  # heme present before

    # Apply just the cofactor pseudoreaction step
    conditions._apply_cofactor_pseudoreaction(
        mutated, conditions.load_condition("anaerobic")
    )
    # heme coefficient should now be zero (or removed)
    assert cofac.metabolites.get(heme, 0) == 0


def test_anaerobic_biomass_delta_adds_fadh2(model):
    """The biomass delta block adds 0.08 FADH2, -0.08 FAD, -0.16 H+ to r_4041."""
    mutated = model.copy()
    bio = mutated.reactions.get_by_id("r_4041")
    fadh2 = mutated.metabolites.get_by_id("s_0689")
    fad = mutated.metabolites.get_by_id("s_0687")
    proton = mutated.metabolites.get_by_id("s_0794")

    before = {
        fadh2.id: bio.metabolites.get(fadh2, 0),
        fad.id: bio.metabolites.get(fad, 0),
        proton.id: bio.metabolites.get(proton, 0),
    }

    conditions._apply_biomass_stoichiometry_delta(
        mutated, conditions.load_condition("anaerobic")
    )

    after = bio.metabolites
    assert after[fadh2] == pytest.approx(before[fadh2.id] + 0.08)
    assert after[fad] == pytest.approx(before[fad.id] - 0.08)
    assert after[proton] == pytest.approx(before[proton.id] - 0.16)
