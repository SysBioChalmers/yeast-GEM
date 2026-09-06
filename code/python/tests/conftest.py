"""Shared pytest fixtures for yeastgem tests."""
from __future__ import annotations

import pytest

from yeastgem import load_yeast_yaml


@pytest.fixture(scope="session")
def model():
    """The yeast-GEM model loaded once per test session."""
    return load_yeast_yaml()


@pytest.fixture
def isolated_paths(model, tmp_path, monkeypatch):
    """Redirect MODEL_PATH, YAML_PATH and the ΔG tsv paths into a tmp
    directory, so write-side tests never touch the real repository files.

    Returns the temp directory (with an empty model/ subdirectory already
    created).
    """
    from yeastgem import io as yio
    from yeastgem import missing_fields as mf

    model_dir = tmp_path / "model"
    model_dir.mkdir()
    monkeypatch.setattr(yio, "MODEL_PATH", model_dir / "yeast-GEM.xml")
    monkeypatch.setattr(yio, "YAML_PATH", model_dir / "yeast-GEM.yml")
    monkeypatch.setattr(yio, "REPO_PATH", tmp_path)

    # ΔG tsv redirection lives in missing_fields. Nothing currently calls
    # save_delta_g/load_delta_g with no explicit path while under this
    # fixture, but this guarantees it never touches the real repo files
    # if something does.
    monkeypatch.setattr(mf, "_MET_TSV", tmp_path / "met.tsv")
    monkeypatch.setattr(mf, "_RXN_TSV", tmp_path / "rxn.tsv")

    return tmp_path
