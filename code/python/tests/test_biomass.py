"""Smoke tests for ``yeastgem.biomass`` against the real yeast-GEM model.

Unit-level coverage of the generic biomass mechanism lives upstream in
``raven_toolbox/tests/test_biomass.py``; here we just confirm that the
yeast config (loaded from ``data/yeastgem/ids.yml``) talks to the
upstream API correctly on the real ~4100-reaction model.
"""
from __future__ import annotations

import pytest

from yeastgem import biomass


def test_yeast_biomass_config_components(model):
    cfg = biomass.yeast_biomass_config()
    names = {c.name for c in cfg.components}
    # Mirrors ids.yml::biomass_components.
    assert names == {
        "protein", "carbohydrate", "RNA", "DNA",
        "lipid_backbone", "ion", "cofactor",
    }


def test_sum_biomass_components_present(model):
    """Every configured component plus the total must be in the output."""
    out = biomass.sum_biomass(model.copy())
    assert set(out) == {
        "protein", "carbohydrate", "RNA", "DNA",
        "lipid_backbone", "ion", "cofactor", "total",
    }
    assert out["total"] > 0.5  # yeast biomass sums close to 1 g/gDW


def test_sum_biomass_total_within_realistic_range(model):
    out = biomass.sum_biomass(model.copy())
    # yeast-GEM's biomass equation sums to ~1.0 by design.
    assert out["total"] == pytest.approx(1.0, abs=0.05)


def test_scale_biomass_lands_on_target(model):
    mutated = model.copy()
    before = biomass.sum_biomass(mutated)["protein"]
    target = before * 0.9
    biomass.scale_biomass(mutated, "protein", target)
    after = biomass.sum_biomass(mutated)["protein"]
    assert after == pytest.approx(target, rel=1e-6)


def test_scale_biomass_with_balance_keeps_total(model):
    mutated = model.copy()
    biomass.scale_biomass(mutated, "protein", 0.5, balance_out="carbohydrate")
    out = biomass.sum_biomass(mutated)
    assert out["total"] == pytest.approx(1.0, rel=1e-4)


def test_set_gam_scales_cofactor_coefficients(model):
    mutated = model.copy()
    bio = mutated.reactions.get_by_id("r_4041")
    atp_met = next(m for m in bio.metabolites if m.name == "ATP")
    before = bio.metabolites[atp_met]
    new_gam = 80
    biomass.set_gam(mutated, new_gam)
    expected_sign = 1 if before > 0 else -1
    assert bio.metabolites[atp_met] == pytest.approx(expected_sign * new_gam)


def test_change_amino_acid_ratio_anaerobic_changes_protein_stoich(model):
    mutated = model.copy()
    rxn = mutated.reactions.get_by_id("r_4047")
    # Pick a tRNA we know is in the AA file (alanine, s_0404).
    sub = mutated.metabolites.get_by_id("s_0404")
    before = rxn.metabolites.get(sub, 0)

    biomass.change_amino_acid_ratio(mutated, aerobic=False)
    after = rxn.metabolites.get(sub, 0)
    # The aerobic and anaerobic columns differ → coefficient should change.
    assert after != pytest.approx(before)


def test_change_amino_acid_ratio_preserves_protein_mass(model):
    """The function should rescale protein content back to its
    pre-switch value via :func:`scale_biomass`."""
    mutated = model.copy()
    before_protein = biomass.sum_biomass(mutated)["protein"]
    biomass.change_amino_acid_ratio(mutated, aerobic=False)
    after_protein = biomass.sum_biomass(mutated)["protein"]
    assert after_protein == pytest.approx(before_protein, rel=1e-4)
