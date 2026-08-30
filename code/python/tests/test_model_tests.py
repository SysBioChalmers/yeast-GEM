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
    result = model_tests.anaerobic_flux_predictions(anaerobic)
    # Fold error is >= 1 by construction (1.0 is exact); yeast9 lands
    # near 1.05 on the median.
    assert 1.0 <= result.median_fold_error <= result.mean_fold_error
    assert 0.0 <= result.fraction_within_two_fold <= 1.0
    assert 0 <= result.n_unpredicted < result.n_measurements
    assert 0.5 <= result.r2 <= 1.0


def test_plot_anaerobic_returns_exchange_metrics(model):
    """plot_anaerobic compares predicted exchange rates with measurements."""
    anaerobic = model.copy()
    conditions.apply(anaerobic, "anaerobic")
    result = model_tests.plot_anaerobic(anaerobic, plot=False)

    assert result.n_measurements == len(result.results)
    assert result.n_measurements > 0
    # Glycerol, ethanol, CO2 and biomass are all produced under anaerobic
    # conditions, and the table reports magnitudes.
    assert (result.results["predicted"] >= -1e-6).all()
    assert 0.0 <= result.fraction_within_error <= 1.0
    assert result.max_relative_error >= result.mean_relative_error
    # Ammonium uptake is coupled to the plasma membrane ATPase; the
    # measured ratio is near 1, so anything wildly off indicates the
    # wrong reactions are being read.
    assert 0.5 < result.ammonium_per_atpase < 2.0


def test_find_duplicated_rxns_returns_list(model, capsys):
    groups = model_tests.find_duplicated_rxns(model.copy())
    assert isinstance(groups, list)
    # yeast-GEM has historically had some duplicate-pair survivors;
    # the function must finish without error regardless of count.
    captured = capsys.readouterr()
    if groups:
        assert "Name:" in captured.out
