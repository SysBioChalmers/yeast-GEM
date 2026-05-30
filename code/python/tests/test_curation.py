"""Tests for ``yeastgem.curation`` against the real yeast-GEM model.

Unit-level coverage of the generic engine lives upstream in
``raven_python/tests/test_curation.py``; here we just verify the
yeast-GEM shim picks up the ``s_`` / ``r_`` prefixes and that real
v8_*/v9_* TSVs apply cleanly.
"""
from __future__ import annotations

import pandas as pd

from yeastgem import curation
from yeastgem.io import REPO_PATH


def test_new_met_uses_s_prefix(model):
    mutated = model.copy()
    df = pd.DataFrame([
        {"metNames": "test_metabolite_phase_6", "comps": "c",
         "formula": "C2H6O", "charge": 0, "inchi": "", "metNotes": ""},
    ])
    result = curation.curate_mets_rxns_genes(mutated, mets_df=df)
    assert len(result.added_metabolites) == 1
    assert result.added_metabolites[0].startswith("s_")


def test_new_rxn_uses_r_prefix(model):
    mutated = model.copy()
    # Use an existing yeast met (s_0794 = H+[c]) to avoid the
    # add-new-met machinery.
    atp = next(m for m in mutated.metabolites if m.name == "ATP" and m.compartment == "c")
    rxns_df = pd.DataFrame([
        {"rxnNames": "phase6 test rxn", "grRules": "", "lb": 0, "ub": 1000,
         "rev": 0, "subSystems": "", "eccodes": "", "rxnNotes": "",
         "rxnReferences": "", "rxnConfidenceScores": ""},
    ])
    coeffs_df = pd.DataFrame([
        {"rxnNames": "phase6 test rxn", "metNames": atp.name, "comps": "c",
         "coefficient": -1.0},
        {"rxnNames": "phase6 test rxn", "metNames": "H+", "comps": "c",
         "coefficient": 1.0},
    ])
    result = curation.curate_mets_rxns_genes(
        mutated, rxns_df=rxns_df, rxns_coeffs_df=coeffs_df,
    )
    assert len(result.added_reactions) == 1
    assert result.added_reactions[0].startswith("r_")


def test_real_curation_tsvs_v8_6_3_volpolyp(model):
    """Apply the v8_6_3 VolPolyP curation files end-to-end. Mostly a
    smoke test: confirm no exception, and that some entities were
    added/updated."""
    mutated = model.copy()
    data_dir = REPO_PATH / "data" / "modelCuration" / "v8_6_3"

    result = curation.curate_mets_rxns_genes_from_tsv(
        mutated,
        mets_tsv=data_dir / "VolPolyPMets.tsv",
        genes_tsv=data_dir / "VolPolyPGenes.tsv",
        rxns_tsv=data_dir / "VolPolyPRxns.tsv",
        rxns_coeffs_tsv=data_dir / "VolPolyPRxnsCoeffs.tsv",
    )
    # We applied a TSV pack — at minimum some entity should land.
    touched = (
        len(result.added_metabolites) + len(result.updated_metabolites)
        + len(result.added_genes) + len(result.updated_genes)
        + len(result.added_reactions) + len(result.updated_reactions)
    )
    assert touched > 0


def test_empty_call_no_op(model):
    mutated = model.copy()
    result = curation.curate_mets_rxns_genes(mutated)
    assert not result
