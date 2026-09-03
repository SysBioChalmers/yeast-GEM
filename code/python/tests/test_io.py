"""Smoke tests for ``yeastgem.io``."""
from __future__ import annotations

import tempfile
import warnings
from pathlib import Path

import cobra

from yeastgem import REPO_PATH, YAML_PATH, load_yeast_yaml, read_yeast_model


def test_repo_path_resolves():
    assert REPO_PATH.is_dir()
    assert (REPO_PATH / "model").is_dir()


def test_yaml_path_exists():
    assert YAML_PATH.exists(), f"yeast-GEM yaml not found at {YAML_PATH}"


def test_model_loads(model):
    assert isinstance(model, cobra.Model)
    # Current model has ~4100 reactions / ~2750 metabolites / ~1140 genes.
    # Use loose lower bounds so the test survives expected growth.
    assert len(model.reactions) > 3000
    assert len(model.metabolites) > 2000
    assert len(model.genes) > 1000


def test_default_load_does_not_apply_bigg_compliance(model):
    assert "x" not in model.compartments  # x = peroxisome under BiGG only


def test_load_yeast_yaml_merges_tsv_annotation(model):
    """load_yeast_yaml must merge in the cross-reference annotation that
    lives in model/{reactions,metabolites,genes}.tsv (yeast-GEM#379), not
    just the yml's own inline fields."""
    annotated = [m for m in model.metabolites if "metanetx.chemical" in m.annotation]
    assert annotated, "no metabolite carries metanetx.chemical -- tsv merge did not run"


def test_bigg_compliant_load_produces_a_valid_model():
    """``make_bigg_compliant=True`` remaps identifiers via a lookup table
    (see ``_make_bigg_compliant``), which is the kind of transform that can
    silently corrupt a model rather than raise -- so check the result is
    still complete and still round-trips through SBML, not just that
    loading it did not crash. Formerly exercised by memote-run.yml, which
    ran the full MEMOTE suite against this conversion on every pull
    request without ever reporting a score anywhere; this is the same
    regression coverage in seconds instead of ~13 minutes.
    """
    model = load_yeast_yaml(make_bigg_compliant=True)
    assert "x" in model.compartments  # x = peroxisome, BiGG-only
    assert len(model.reactions) > 3000
    assert len(model.metabolites) > 2000

    with tempfile.TemporaryDirectory() as tmp_dir:
        sbml_path = str(Path(tmp_dir) / "model.xml")
        cobra.io.write_sbml_model(model, sbml_path)
        reloaded = cobra.io.read_sbml_model(sbml_path)
    assert len(reloaded.reactions) == len(model.reactions)
    assert len(reloaded.metabolites) == len(model.metabolites)


def test_load_yeast_yaml_returns_independent_instances():
    a = load_yeast_yaml()
    b = load_yeast_yaml()
    assert a is not b
    assert len(a.reactions) == len(b.reactions)


def test_read_yeast_model_deprecation_warning():
    """The shim must warn and still return a usable, fully-annotated
    model (forwarding to load_yeast_yaml)."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model = read_yeast_model()
    assert any(
        issubclass(w.category, DeprecationWarning)
        and "load_yeast_yaml" in str(w.message)
        for w in caught
    )
    assert len(model.reactions) > 3000
