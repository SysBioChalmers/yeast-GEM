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
    """Redirect MODEL_PATH, YAML_PATH and the ΔG CSV paths into a tmp
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

    # ΔG CSV redirection lives in missing_fields.
    monkeypatch.setattr(mf, "_MET_CSV", tmp_path / "met.csv")
    monkeypatch.setattr(mf, "_RXN_CSV", tmp_path / "rxn.csv")

    return tmp_path
