"""Smoke tests for ``yeastgem.io``."""
from __future__ import annotations

import cobra

from yeastgem import MODEL_PATH, REPO_PATH, read_yeast_model


def test_repo_path_resolves():
    assert REPO_PATH.is_dir()
    assert (REPO_PATH / "model").is_dir()


def test_model_path_exists():
    assert MODEL_PATH.exists(), f"yeast-GEM SBML not found at {MODEL_PATH}"


def test_model_loads(model):
    assert isinstance(model, cobra.Model)
    # Current model has ~4100 reactions / ~2750 metabolites / ~1140 genes.
    # Use loose lower bounds so the test survives expected growth.
    assert len(model.reactions) > 3000
    assert len(model.metabolites) > 2000
    assert len(model.genes) > 1000


def test_default_load_does_not_apply_bigg_compliance(model):
    assert "x" not in model.compartments  # x = peroxisome under BiGG only


def test_read_yeast_model_returns_independent_instances():
    a = read_yeast_model()
    b = read_yeast_model()
    assert a is not b
    assert len(a.reactions) == len(b.reactions)
