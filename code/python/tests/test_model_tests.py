"""Tests for yeastgem.model_tests against the real yeast-GEM model.

Slow tests — each exercises FBA on the full 4000+ reaction model.
Tolerances are deliberately loose; the strict pass/fail thresholds
live in the lock-step verification driver (see PORTING_PLAN.md and
tests/reference/).
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from yeastgem import conditions, model_tests


def test_growth_returns_high_r2(model):
    """Yeast8/9 chemostat R² is ~0.99. Anything ≥ 0.9 means the
    chemostat sweep is configured correctly."""
    r2 = model_tests.growth(model.copy())
    assert 0.9 <= r2 <= 1.0


def test_growth_returns_float(model):
    r2 = model_tests.growth(model.copy())
    assert isinstance(r2, float)


def test_essential_genes_returns_reasonable_accuracy(model):
    """Yeast9 typically lands around 0.84 accuracy with 700+ verified
    genes covered. Anything above 0.7 means cobra single_gene_deletion
    ran and the Stanford reference lists loaded."""
    result = model_tests.essential_genes(model.copy())
    assert 0.7 <= result.accuracy <= 1.0
    assert len(result.tp) + len(result.tn) + len(result.fp) + len(result.fn) > 500


def test_essential_genes_sensitivity_specificity(model):
    """Both rates must be in [0, 100] (percent)."""
    result = model_tests.essential_genes(model.copy())
    assert 0 <= result.sensitivity <= 100
    assert 0 <= result.specificity <= 100


def test_anaerobic_flux_predictions(model):
    """Apply the anaerobic condition first; the function expects it
    pre-applied (mirrors the legacy MATLAB calling convention)."""
    anaerobic = model.copy()
    conditions.apply(anaerobic, "anaerobic")
    r2, mre = model_tests.anaerobic_flux_predictions(anaerobic)
    assert 0.5 <= r2 <= 1.0  # yeast9 typically lands ~0.95
    assert 0 <= mre  # MRE is non-negative


def test_plot_anaerobic_returns_predictions(model):
    """plot_anaerobic also returns the (gly, eth, CO2, biomass) vector."""
    anaerobic = model.copy()
    conditions.apply(anaerobic, "anaerobic")
    sim = model_tests.plot_anaerobic(anaerobic)
    assert sim.shape == (4,)
    # Glycerol, ethanol, CO2 should all be non-negative under anaerobic;
    # biomass too (growth rate).
    assert (sim >= -1e-6).all()


def test_find_duplicated_rxns_returns_list(model, capsys):
    groups = model_tests.find_duplicated_rxns(model.copy())
    assert isinstance(groups, list)
    # yeast-GEM has historically had some duplicate-pair survivors;
    # the function must finish without error regardless of count.
    captured = capsys.readouterr()
    if groups:
        assert "Name:" in captured.out
