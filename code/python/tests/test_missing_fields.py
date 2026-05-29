"""Tests for ``yeastgem.missing_fields`` (SBO + ΔG persistence)."""
from __future__ import annotations

import math

import pandas as pd

from yeastgem import add_sbo_terms, load_delta_g, save_delta_g
from yeastgem.missing_fields import (
    _DELTA_G_NOTE_KEY,
    _transport_reaction_ids,
)

# --- add_sbo_terms ---------------------------------------------------

def test_add_sbo_terms_assigns_every_met(model):
    mutated = model.copy()
    add_sbo_terms(mutated)
    for met in mutated.metabolites:
        assert met.annotation.get("sbo"), f"{met.id} has no SBO term"


def test_add_sbo_terms_assigns_every_rxn(model):
    mutated = model.copy()
    add_sbo_terms(mutated)
    for rxn in mutated.reactions:
        assert rxn.annotation.get("sbo"), f"{rxn.id} has no SBO term"


def test_biomass_pseudo_metabolites_get_biomass_sbo(model):
    mutated = model.copy()
    add_sbo_terms(mutated)
    biomass_mets = [m for m in mutated.metabolites if m.name == "biomass"]
    assert biomass_mets, "no metabolite named 'biomass' found in test fixture"
    for met in biomass_mets:
        assert met.annotation["sbo"] == "SBO:0000649"


def test_simple_chemicals_get_simple_chemical_sbo(model):
    mutated = model.copy()
    add_sbo_terms(mutated)
    atp_mets = [m for m in mutated.metabolites if m.name == "ATP"]
    assert atp_mets
    for met in atp_mets:
        assert met.annotation["sbo"] == "SBO:0000247"


def test_exchange_reactions_get_exchange_sbo(model):
    mutated = model.copy()
    add_sbo_terms(mutated)
    # Exchange reactions in cobra: single met in extracellular compartment.
    exchanges = list(mutated.exchanges)
    assert exchanges, "no exchanges in test fixture"
    extracellular_exchanges = [
        r for r in exchanges
        if len(r.metabolites) == 1 and next(iter(r.metabolites)).compartment == "e"
    ]
    for rxn in extracellular_exchanges:
        assert rxn.annotation["sbo"] == "SBO:0000627", (
            f"{rxn.id} got {rxn.annotation['sbo']}"
        )


def test_fill_semantic_preserves_existing(model):
    mutated = model.copy()
    met = mutated.metabolites[0]
    met.annotation["sbo"] = "SBO:0009999"  # arbitrary pre-existing value
    add_sbo_terms(mutated)
    assert met.annotation["sbo"] == "SBO:0009999"


def test_transport_reaction_detection(model):
    """A reaction with the same metabolite name in two compartments
    should be flagged as transport."""
    transport_ids = _transport_reaction_ids(model)
    assert transport_ids, "expected at least some transport reactions"
    # Sanity: an exchange reaction is NOT a transport reaction.
    exchange_ids = {r.id for r in model.exchanges}
    assert not (transport_ids & exchange_ids)


# --- ΔG load / save round-trip --------------------------------------

def test_load_delta_g_populates_notes(model):
    mutated = model.copy()
    load_delta_g(mutated)
    stamped = sum(
        1 for m in mutated.metabolites if _DELTA_G_NOTE_KEY in m.notes
    )
    assert stamped > 0, "expected load_delta_g to stamp at least some metabolites"


def test_save_delta_g_round_trip(model, tmp_path):
    mutated = model.copy()
    load_delta_g(mutated)

    met_csv = tmp_path / "met.csv"
    rxn_csv = tmp_path / "rxn.csv"
    save_delta_g(mutated, met_csv=met_csv, rxn_csv=rxn_csv)

    assert met_csv.exists() and rxn_csv.exists()
    met_df = pd.read_csv(met_csv)
    assert list(met_df.columns) == ["Var1", "Var2"]
    assert len(met_df) == len(mutated.metabolites)
    assert list(met_df["Var1"]) == [m.id for m in mutated.metabolites]


def test_save_then_load_preserves_values(model, tmp_path):
    """Save ΔG, reload into a fresh model, and verify the notes survive."""
    seed = model.copy()
    load_delta_g(seed)
    met_csv = tmp_path / "met.csv"
    rxn_csv = tmp_path / "rxn.csv"
    save_delta_g(seed, met_csv=met_csv, rxn_csv=rxn_csv)

    fresh = model.copy()
    load_delta_g(fresh, met_csv=met_csv, rxn_csv=rxn_csv)

    for original, reloaded in zip(seed.metabolites, fresh.metabolites, strict=True):
        assert original.notes.get(_DELTA_G_NOTE_KEY) == \
            reloaded.notes.get(_DELTA_G_NOTE_KEY)


def test_save_delta_g_emits_nan_for_missing_notes(model, tmp_path):
    """Metabolites without a ΔG note appear in the CSV as NaN, preserving
    one-row-per-entity ordering (mirrors MATLAB's array2table behaviour).

    Note: the committed SBML already carries ΔG values in metabolite
    notes from MATLAB's release pipeline, so we must explicitly clear
    them before asserting the "missing → NaN" behaviour.
    """
    fresh = model.copy()
    for met in fresh.metabolites:
        met.notes.pop(_DELTA_G_NOTE_KEY, None)
    for rxn in fresh.reactions:
        rxn.notes.pop(_DELTA_G_NOTE_KEY, None)

    met_csv = tmp_path / "met.csv"
    rxn_csv = tmp_path / "rxn.csv"
    save_delta_g(fresh, met_csv=met_csv, rxn_csv=rxn_csv)
    df = pd.read_csv(met_csv)
    assert len(df) == len(fresh.metabolites)
    assert df["Var2"].apply(lambda v: math.isnan(v)).all()
