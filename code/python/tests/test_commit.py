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

    Copies the canonical README so the regex rewrite has something
    to chew on. Returns the temp directory.
    """
    from yeastgem import missing_fields as mf

    model_dir = tmp_path / "model"
    model_dir.mkdir()
    model_path = model_dir / "yeast-GEM.xml"
    monkeypatch.setattr(yio, "MODEL_PATH", model_path)
    monkeypatch.setattr(yio, "REPO_PATH", tmp_path)

    # README seed: the regex looks for a yeast-GEM stats row.
    (tmp_path / "README.md").write_text(
        "Header\n"
        "| Taxonomy | Latest update | Version | Reactions | Metabolites | Genes |\n"
        "|:-------|:--------------|:------|:------|:----------|:-----|\n"
        "| _Saccharomyces cerevisiae_ | 01-Jan-2000 | develop | 1 | 1 | 1 |\n"
        "Footer\n",
        encoding="utf-8",
    )

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


def test_commit_yeast_model_updates_readme(model, isolated_paths):
    commit_yeast_model(model.copy())
    text = (isolated_paths / "README.md").read_text(encoding="utf-8")
    # Old stub row was 1/1/1; the rewrite plugs in real model sizes.
    assert "| 1 | 1 | 1 |" not in text
    assert "| _Saccharomyces cerevisiae_" in text
    # Today's date should appear (just check the year)
    from datetime import datetime
    assert datetime.now().strftime("%Y") in text


def test_commit_yeast_model_skip_readme(model, isolated_paths):
    commit_yeast_model(model.copy(), update_readme=False)
    text = (isolated_paths / "README.md").read_text(encoding="utf-8")
    assert "01-Jan-2000" in text  # stub preserved


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


def test_commit_anaerobic_check_warns(model, isolated_paths):
    """Anaerobic growth check is deferred — must emit a warning by default."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        commit_yeast_model(model.copy())
    assert any("Anaerobic growth check is deferred" in str(w.message) for w in caught)


def test_commit_anaerobic_strict_raises(model, isolated_paths):
    """With allow_no_growth=False the anaerobic gap must raise."""
    with pytest.raises(NotImplementedError, match="Anaerobic growth check"):
        commit_yeast_model(model.copy(), allow_no_growth=False)


def test_commit_returns_model(model, isolated_paths):
    mutated = model.copy()
    returned = commit_yeast_model(mutated)
    assert returned is mutated
    assert isinstance(returned, cobra.Model)
