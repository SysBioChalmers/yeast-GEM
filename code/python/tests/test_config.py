"""Tests for ``yeastgem.config`` (canonical yeast IDs)."""
from __future__ import annotations

from yeastgem import YeastIDs, load_ids


def test_load_ids_returns_dataclass():
    ids = load_ids()
    assert isinstance(ids, YeastIDs)


def test_load_ids_has_expected_keys():
    ids = load_ids()
    assert ids.biomass_rxn == "r_4041"
    assert ids.protein_rxn == "r_4047"
    assert ids.cofactor_rxn == "r_4598"
    assert ids.proton_met == "s_0794"


def test_load_ids_pseudoreaction_map_complete():
    ids = load_ids()
    expected = {
        "biomass", "protein", "carbohydrate",
        "lipid_backbone", "lipid_chain",
        "RNA", "DNA", "ion", "cofactor",
    }
    assert set(ids.pseudoreaction_names) == expected
    assert ids.pseudoreaction_names["biomass"] == "biomass pseudoreaction"


def test_load_ids_gam_cofactors_match_legacy():
    ids = load_ids()
    # Mirrors the hardcoded list in legacy changeGAM.m
    assert ids.gam_cofactors == ["ATP", "ADP", "H2O", "H+", "phosphate"]


def test_ids_exist_in_model(model):
    ids = load_ids()
    assert ids.biomass_rxn in {r.id for r in model.reactions}
    assert ids.protein_rxn in {r.id for r in model.reactions}
    assert ids.cofactor_rxn in {r.id for r in model.reactions}
    assert ids.proton_met in {m.id for m in model.metabolites}
