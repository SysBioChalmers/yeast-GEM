"""Tests for ``commit_yeast_model`` and the ``write_yeast_model`` shim."""
from __future__ import annotations

import warnings

import cobra
import pytest

from yeastgem import commit_yeast_model, write_yeast_model
from yeastgem import io as yio


def test_write_yeast_model_deprecation_warning(model, isolated_paths):
    """The shim must warn and still write a usable model."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        write_yeast_model(model.copy())
    assert any(
        issubclass(w.category, DeprecationWarning)
        and "commit_yeast_model" in str(w.message)
        for w in caught
    )
    assert yio.MODEL_PATH.exists()


def test_commit_yeast_model_writes_sbml(model, isolated_paths):
    commit_yeast_model(model.copy())
    assert yio.MODEL_PATH.exists()
    assert yio.MODEL_PATH.stat().st_size > 1_000_000  # ~MB-scale SBML


def test_commit_yeast_model_writes_txt_by_default(model, isolated_paths):
    commit_yeast_model(model.copy())
    assert (yio.MODEL_PATH.parent / "yeast-GEM.txt").exists()


def test_commit_yeast_model_formats_xml_only(model, isolated_paths):
    commit_yeast_model(model.copy(), formats=("xml",))
    assert yio.MODEL_PATH.exists()
    assert not (yio.MODEL_PATH.parent / "yeast-GEM.txt").exists()


def test_commit_yeast_model_does_not_write_yaml(model, isolated_paths):
    """commit_yeast_model must not touch model/yeast-GEM.yml or the
    annotation tsvs -- that is exclusively save_yeast_yaml's job."""
    commit_yeast_model(model.copy())
    assert not yio.YAML_PATH.exists()


def test_commit_yeast_model_rejects_unsupported_formats(model, isolated_paths):
    """'xlsx'/'mat' need RAVEN; only 'xml'/'txt' are available in Python."""
    with pytest.raises(ValueError, match="xlsx"):
        commit_yeast_model(model.copy(), formats=("xlsx",))


def test_commit_yeast_model_never_ships_delta_g(model, isolated_paths):
    """ΔG is an estimated, not curator-verified value: it must never
    appear in the exported model file, even if the caller's model
    happens to carry it, and commit_yeast_model must never touch the
    ΔG tsvs at all (call load_delta_g/save_delta_g explicitly for that)."""
    from yeastgem import missing_fields as mf

    mutated = model.copy()
    # Stamp a deltaG note directly rather than going through load_delta_g:
    # this test only needs "the model happens to carry deltaG", not a
    # real load, and load_delta_g's tsv paths default to
    # isolated_paths' redirected (deliberately nonexistent) tmp files.
    met = mutated.metabolites[0]
    met.notes = {**met.notes, "deltaG": "-343.18"}

    commit_yeast_model(mutated)

    assert not mf._MET_TSV.exists() and not mf._RXN_TSV.exists()

    reloaded = cobra.io.read_sbml_model(str(yio.MODEL_PATH))
    assert not any("deltaG" in m.notes for m in reloaded.metabolites)
    # commit_yeast_model must not have stripped the caller's own model.
    assert any("deltaG" in m.notes for m in mutated.metabolites)


def test_commit_applies_canonical_state(model, isolated_paths):
    """After commit, the model must have minimal_Y6 bounds + SBO annotations."""
    mutated = model.copy()
    commit_yeast_model(mutated)
    # Bicarbonate exchange should be blocked (minimal_Y6).
    bicarb = mutated.reactions.get_by_id("r_1663")
    assert bicarb.lower_bound == 0 and bicarb.upper_bound == 0
    # Every reaction has an SBO annotation.
    for rxn in mutated.reactions:
        assert rxn.annotation.get("sbo")


def test_commit_runs_anaerobic_growth_check(model, isolated_paths):
    """Phase 4 turned the anaerobic check on. The deferred-warning
    message must no longer appear (and the pipeline must finish)."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        commit_yeast_model(model.copy())
    assert not any("deferred to phase" in str(w.message) for w in caught)


def test_commit_anaerobic_strict_succeeds(model, isolated_paths):
    """With allow_no_growth=False the anaerobic check now runs end-to-end."""
    commit_yeast_model(model.copy(), allow_no_growth=False)
    # Pipeline returns normally; bombing out would have raised.


def test_commit_returns_model(model, isolated_paths):
    mutated = model.copy()
    returned = commit_yeast_model(mutated)
    assert returned is mutated
    assert isinstance(returned, cobra.Model)
