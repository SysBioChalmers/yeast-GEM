"""Tests for the level-1 model comparator."""
from __future__ import annotations

import pytest

from yeastgem import compare_models, read_yeast_model
from yeastgem.compare import _normalise_gpr


def test_model_equal_to_itself(model):
    report = compare_models(model, model)
    assert report.equal, report


def test_independent_loads_are_equal():
    a = read_yeast_model()
    b = read_yeast_model()
    report = compare_models(a, b)
    assert report.equal, report


def test_dropped_reaction_is_detected(model):
    mutated = model.copy()
    rxn = mutated.reactions[0]
    mutated.remove_reactions([rxn])
    report = compare_models(model, mutated)
    assert not report.equal
    assert any("reactions only in A" in d for d in report.differences)


def test_bound_difference_is_detected(model):
    mutated = model.copy()
    rxn = next(r for r in mutated.reactions if r.lower_bound != -42.0)
    rxn.lower_bound = -42.0
    report = compare_models(model, mutated)
    assert not report.equal
    assert any("bounds" in d for d in report.differences)


def test_stoichiometry_within_tolerance_is_equal(model):
    mutated = model.copy()
    rxn = mutated.reactions[0]
    met = next(iter(rxn.metabolites))
    rxn.add_metabolites({met: 1e-12}, combine=True)  # nudge within tol
    report = compare_models(model, mutated, stoichiometry_tol=1e-9)
    assert report.equal, report


def test_stoichiometry_outside_tolerance_is_not_equal(model):
    mutated = model.copy()
    rxn = mutated.reactions[0]
    met = next(iter(rxn.metabolites))
    rxn.add_metabolites({met: 999.0}, combine=False)
    report = compare_models(model, mutated, stoichiometry_tol=1e-9)
    assert not report.equal
    assert any("coef[" in d for d in report.differences)


@pytest.mark.parametrize(
    "a, b",
    [
        ("A and B", "a AND b"),
        ("A  and   B", "A and B"),
        ("(A or B) and C", "(a OR b) AND c"),
    ],
)
def test_gpr_normalisation_is_case_and_whitespace_insensitive(a, b):
    assert _normalise_gpr(a) == _normalise_gpr(b)


def test_report_str_lists_differences(model):
    mutated = model.copy()
    mutated.remove_reactions([mutated.reactions[0]])
    report = compare_models(model, mutated)
    text = str(report)
    assert "Models differ" in text
    assert "reactions only in A" in text
