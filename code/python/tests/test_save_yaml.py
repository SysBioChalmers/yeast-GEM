"""Tests for ``save_yeast_yaml``."""
from __future__ import annotations

from raven_toolbox.io import read_yaml_model

from yeastgem import io as yio
from yeastgem import save_yeast_yaml


def test_save_yeast_yaml_writes_yml(model, isolated_paths):
    save_yeast_yaml(model.copy())
    assert yio.YAML_PATH.exists()


def test_save_yeast_yaml_does_not_write_binaries(model, isolated_paths):
    """save_yeast_yaml must never touch SBML/xml/txt/xlsx/mat -- that is
    exclusively commit_yeast_model's job."""
    save_yeast_yaml(model.copy())
    assert not yio.MODEL_PATH.exists()


def test_save_yeast_yaml_writes_deltag_csvs(model, isolated_paths):
    from yeastgem import missing_fields as mf

    save_yeast_yaml(model.copy())
    assert mf._MET_CSV.exists() and mf._RXN_CSV.exists()


def test_save_yeast_yaml_applies_canonical_state(model, isolated_paths):
    """After save, the model must have minimal_Y6 bounds + SBO annotations."""
    mutated = model.copy()
    save_yeast_yaml(mutated)
    # Bicarbonate exchange should be blocked (minimal_Y6).
    bicarb = mutated.reactions.get_by_id("r_1663")
    assert bicarb.lower_bound == 0 and bicarb.upper_bound == 0
    for rxn in mutated.reactions:
        assert rxn.annotation.get("sbo")


def test_save_yeast_yaml_strips_tsv_annotation_from_written_yml(model, isolated_paths):
    """model/yeast-GEM.yml must not re-embed the cross-reference
    annotation that lives in the tsvs (yeast-GEM#379) -- but save_yeast_yaml
    must leave the caller's own in-memory model untouched, since it only
    strips a copy before writing."""
    mutated = model.copy()
    assert any("metanetx.chemical" in m.annotation for m in mutated.metabolites)

    save_yeast_yaml(mutated)

    assert any("metanetx.chemical" in m.annotation for m in mutated.metabolites), \
        "save_yeast_yaml must not strip the caller's own model"

    written = read_yaml_model(str(yio.YAML_PATH))
    assert not any("metanetx.chemical" in m.annotation for m in written.metabolites), \
        "the written yml must not re-embed tsv-sourced annotation"


def test_save_yeast_yaml_returns_model(model, isolated_paths):
    mutated = model.copy()
    returned = save_yeast_yaml(mutated)
    assert returned is mutated


def test_save_yeast_yaml_derive_tsvs_off_by_default(model, isolated_paths):
    """derive_tsvs defaults to False: a curator's hand-edited tsv must
    never be silently overwritten."""
    reactions_tsv = yio.YAML_PATH.parent / "reactions.tsv"
    header = "id\tname\tbigg.reaction\tec-code\tkegg.pathway\tkegg.reaction\tmetanetx.reaction"
    sentinel = f"{header}\nSENTINEL\n"
    reactions_tsv.write_text(sentinel, encoding="utf-8")

    save_yeast_yaml(model.copy())

    assert reactions_tsv.read_text(encoding="utf-8") == sentinel


def test_save_yeast_yaml_derive_tsvs_writes_programmatic_annotation(model, isolated_paths):
    """A cross-reference added programmatically (not via a tsv edit) must
    survive a derive_tsvs=True save -- the gap deriveAnnotationTsvs closes."""
    mutated = model.copy()
    met = mutated.metabolites[0]
    # Reassign a new dict rather than mutating met.annotation in place:
    # cobra.Model.copy() does not deep-copy per-entity annotation dicts,
    # so an in-place edit here would also alter the session-scoped model
    # fixture's own metabolite, leaking into later tests.
    met.annotation = {**met.annotation, "chebi": "CHEBI:99999999"}

    save_yeast_yaml(mutated, derive_tsvs=True)

    written = (yio.YAML_PATH.parent / "metabolites.tsv").read_text(encoding="utf-8")
    row = next(line for line in written.splitlines() if line.startswith(f"{met.id}\t"))
    assert "CHEBI:99999999" in row
