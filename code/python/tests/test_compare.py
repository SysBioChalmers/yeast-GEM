"""Smoke tests of the cross-language comparator against the real yeast-GEM model.

Unit-level coverage of ``diff_models`` lives upstream in
``raven_toolbox/tests/test_comparison_diff.py``; here we just verify the
yeast-GEM shim re-exports correctly and that the comparator behaves
sensibly on the real 4000+ reaction model.
"""
from __future__ import annotations

from yeastgem import ComparisonReport, compare_models, read_yeast_model


def test_real_model_equal_to_itself(model):
    report = compare_models(model, model)
    assert isinstance(report, ComparisonReport)
    assert report.equal, report


def test_real_model_independent_loads_are_equal():
    a = read_yeast_model()
    b = read_yeast_model()
    assert compare_models(a, b).equal


def test_real_model_dropped_reaction_is_detected(model):
    mutated = model.copy()
    mutated.remove_reactions([mutated.reactions[0]])
    report = compare_models(model, mutated)
    assert not report.equal
    assert any("reactions only in A" in d for d in report.differences)
