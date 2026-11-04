"""Tests for the pull-request report and the pipeline steps around it."""
from __future__ import annotations

import io
import urllib.error

import pytest

from tests.ci import build_qc_report as bq
from tests.ci import pipeline

BASE = {"growth": 0.08, "orphan_metabolites": 3, "duplicate_keys": 0, "empty_reactions": 0,
        "mass_imbalanced": 8, "n_reactions": 100}


def _report(current, base=BASE, status=None, running=False):
    return bq.collect(current, base, "develop", running, status=status)


def _condensed(report, current, base=BASE, running=False):
    return bq.render_condensed(report, current, base, "develop", "https://x/u", running,
                               "https://x/d", "https://x/d/model_qc_summary.md")


def test_regression_is_listed_and_counted():
    current = {**BASE, "mass_imbalanced": 10}
    report = _report(current)
    text = _condensed(report, current)
    assert report["regressions"] == 1
    assert "**Regressions**" in text
    assert "Mass-imbalanced reactions" in text.split("**Regressions**")[1].split("**")[0]
    assert "1 regression(s)" in text


def test_improvement_and_pre_existing_are_listed_apart():
    current = {**BASE, "orphan_metabolites": 1}
    text = _condensed(_report(current), current)
    improved = text.split("**Improved**")[1].split("**")[0]
    assert "Unused metabolites" in improved
    assert "Mass-imbalanced reactions" in text.split("**Pre-existing**")[1]


def test_a_failed_status_check_blocks_the_merge():
    report = _report(BASE, status={"roundtrip_yaml": "pass", "roundtrip_sbml": "fail",
                                   "yamllint": "pass"})
    assert report["fatal_keys"] == ["roundtrip_sbml"]
    assert "the SBML round-trip failed" in bq._verdict(report, BASE, "develop")


def test_a_check_without_a_value_is_unknown_not_clean():
    current = {k: v for k, v in BASE.items() if k != "mass_imbalanced"}
    report = _report(current)
    assert report["skipped"] >= 1
    assert "did not run" in bq._verdict(report, BASE, "develop")


def test_running_report_marks_every_row_as_running():
    report = _report({}, running=True)
    text = _condensed(report, {}, running=True)
    assert report["pending"] > 0
    assert "running" in text and ":white_check_mark:" not in text.split("Section")[1]


def test_full_report_has_every_row_and_no_run_link():
    current = {**BASE, "mass_imbalanced": 10}
    report = _report(current)
    full = bq.render_full(report, current, BASE, "develop", "https://x/u", False, "",
                          "https://x/d")
    assert "[Mass-imbalanced reactions](https://x/u/README.md#mass-imbalanced-reactions)" in full
    assert "https://x/d/qc_findings.md#mass-imbalanced-reactions" in full
    assert "Full workflow run" not in full


def _write(path, text):
    path.write_text(text, encoding="utf-8")
    return path


def test_gate_fails_on_a_duplicate_key(tmp_path):
    metrics = _write(tmp_path / "m.tsv", "growth\t0.1\nduplicate_keys\t1\n")
    assert pipeline.main(["gate", "--metrics", str(metrics)]) == 1


def test_gate_reads_the_file_status(tmp_path):
    metrics = _write(tmp_path / "m.tsv", "growth\t0.1\n")
    good = _write(tmp_path / "good.tsv",
                  "check\tresult\nroundtrip_sbml\tpass\nroundtrip_yaml\tpass\nyamllint\tpass\n")
    bad = _write(tmp_path / "bad.tsv",
                 "check\tresult\nroundtrip_sbml\tpass\nroundtrip_yaml\tfail\nyamllint\tpass\n")
    assert pipeline.main(["gate", "--metrics", str(metrics), "--status", str(good)]) == 0
    assert pipeline.main(["gate", "--metrics", str(metrics), "--status", str(bad)]) == 1
    missing = tmp_path / "none.tsv"
    assert pipeline.main(["gate", "--metrics", str(metrics), "--status", str(missing)]) == 1


def test_stamp_writes_the_pull_request_number(tmp_path):
    readme = _write(tmp_path / "README.md",
                    "| a | <!-- last:qc -->**PR #1** |\n| b | <!-- last:memote -->**PR #2** |\n")
    assert pipeline.main(["stamp", "--readme", str(readme), "--pr", "431",
                          "--group", "qc"]) == 0
    assert readme.read_text(encoding="utf-8") == (
        "| a | <!-- last:qc -->**PR #431** |\n| b | <!-- last:memote -->**PR #2** |\n")


def test_request_retries_a_transient_error(monkeypatch):
    calls = []

    class Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    def fake_urlopen(request):
        calls.append(request.full_url)
        if len(calls) < 3:
            raise urllib.error.HTTPError(request.full_url, 503, "Unavailable", {}, None)
        return Response(b'{"id": 7}')

    monkeypatch.setattr(pipeline.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(pipeline.time, "sleep", lambda seconds: None)
    assert pipeline._request("GET", "https://api.github.com/x", "token") == {"id": 7}
    assert len(calls) == 3


def test_request_does_not_retry_a_client_error(monkeypatch):
    calls = []

    def fake_urlopen(request):
        calls.append(1)
        raise urllib.error.HTTPError(request.full_url, 404, "Not Found", {}, None)

    monkeypatch.setattr(pipeline.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(pipeline.time, "sleep", lambda seconds: None)
    with pytest.raises(urllib.error.HTTPError):
        pipeline._request("GET", "https://api.github.com/x", "token")
    assert len(calls) == 1


# --- gene essentiality ------------------------------------------------

NAMES = {"YEL047C": "FRD1", "YJR051W": "OSM1", "YAL001C": "TFC3"}


def test_gene_changes_are_classified_by_effect():
    base = {"YEL047C": "FP", "YJR051W": "TP", "YAL001C": "TN", "YAL002W": "FN", "YAL003W": "TP"}
    now = {"YEL047C": "TN", "YJR051W": "FN", "YAL001C": "TN", "YAL002W": "FP", "YAL004W": "TP"}
    summary = bq.essential_summary(now, base, NAMES)
    changed = {row[0]: row[4] for row in summary["changed"]}
    assert changed == {"YEL047C": "better", "YJR051W": "worse", "YAL002W": "changed",
                       "YAL003W": "removed", "YAL004W": "new"}
    assert summary["false_non_essential"] == [("YAL002W", "")]
    assert summary["false_essential"] == [("YJR051W", "OSM1")]


def test_without_calls_or_baseline_nothing_is_claimed():
    assert bq.essential_summary({}, {"YEL047C": "TP"}, NAMES) is None
    summary = bq.essential_summary({"YEL047C": "FP"}, {}, NAMES)
    assert summary["changed"] == [] and not summary["has_base"]
    assert "no baseline" in bq._essential_cell(summary, "develop", False)


def test_the_comment_names_the_changed_genes_and_counts_the_wrong_calls():
    base = {"YEL047C": "FP", "YJR051W": "TP"}
    now = {"YEL047C": "TN", "YJR051W": "FN"}
    essential = bq.essential_summary(now, base, NAMES)
    report = bq.collect(BASE, BASE, "develop", False, essential=essential)
    text = bq.render_condensed(report, BASE, BASE, "develop", "https://x/u", False,
                               "https://x/d", "https://x/d/model_qc_summary.md")
    assert "| Gene essentiality | :warning: 1 better, 1 worse vs `develop`; " \
           "0 false non-essential, 1 false essential |" in text
    assert "- `YEL047C` (FRD1): false non-essential → true essential (better)" in text
    assert "- `YJR051W` (OSM1): true non-essential → false essential (worse)" in text


def test_the_comment_lists_at_most_twenty_changed_genes():
    base = {f"G{i:03d}": "FP" for i in range(30)}
    now = {gene: "TN" for gene in base}
    essential = bq.essential_summary(now, base, {})
    report = bq.collect(BASE, BASE, "develop", False, essential=essential)
    text = bq.render_condensed(report, BASE, BASE, "develop", "https://x/u", False,
                               "https://x/d", "https://x/d/model_qc_summary.md")
    assert text.count("(better)") == 20
    assert "and 10 more, see the [full report](https://x/d/model_qc_summary.md#gene-essentiality)" \
        in text


def test_the_full_report_lists_every_change_and_every_wrong_call():
    base = {"YEL047C": "FP"}
    now = {"YEL047C": "TN", "YAL001C": "FN"}
    essential = bq.essential_summary(now, base, NAMES)
    report = bq.collect(BASE, BASE, "develop", False, essential=essential)
    full = bq.render_full(report, BASE, BASE, "develop", "https://x/u", False, "",
                          "https://x/d")
    assert "### Gene essentiality" in full
    assert "| `YEL047C` | FRD1 | false non-essential | true essential | better |" in full
    assert "#### False essential genes (1)" in full and "| `YAL001C` | TFC3 |" in full


def test_gene_calls_and_names_are_read_from_their_files(tmp_path):
    (tmp_path / "essentialGenes.tsv").write_text("gene\tclassification\nYEL047C\tTP\nYAL001C\tFN\n")
    (tmp_path / "names.tsv").write_text("Systematic_name\tStandard_name\nYEL047C\tFRD1\nQ0045\t\n")
    assert bq.read_calls(tmp_path) == {"YEL047C": "TP", "YAL001C": "FN"}
    assert bq.read_calls(tmp_path / "missing") == {}
    assert bq.read_gene_names(tmp_path / "names.tsv") == {"YEL047C": "FRD1"}
