"""Tests for ``commit_yeast_model`` and the ``write_yeast_model`` shim."""
from __future__ import annotations

import warnings

import cobra
import pytest

from yeastgem import commit_yeast_model, write_yeast_model
from yeastgem import io as yio


@pytest.fixture
def isolated_paths(model, tmp_path, monkeypatch):
    """Redirect MODEL_PATH and ΔG CSV paths into a tmp directory.

    Returns the temp directory.
    """
    from yeastgem import missing_fields as mf

    model_dir = tmp_path / "model"
    model_dir.mkdir()
    model_path = model_dir / "yeast-GEM.xml"
    monkeypatch.setattr(yio, "MODEL_PATH", model_path)
    monkeypatch.setattr(yio, "REPO_PATH", tmp_path)

    # ΔG CSV redirection lives in missing_fields.
    monkeypatch.setattr(mf, "_MET_CSV", tmp_path / "met.csv")
    monkeypatch.setattr(mf, "_RXN_CSV", tmp_path / "rxn.csv")

    return tmp_path


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


def test_commit_yeast_model_writes_deltag_csvs(model, isolated_paths):
    """commit pipeline must persist ΔG CSVs to the configured paths."""
    from yeastgem import missing_fields as mf

    commit_yeast_model(model.copy())
    assert mf._MET_CSV.exists() and mf._RXN_CSV.exists()


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
