"""Run the MEMOTE test suite in-process and record the total score.

MEMOTE only computes a score when it builds a report; the raw ``memote run``
result holds per-test results but no score. So this runs MEMOTE through its
Python API and then builds the snapshot report to obtain the score, rather
than shelling out to ``memote run`` (which would leave us without one). One
run per pull request rather than a core-subset/full-suite split -- yeast-GEM
does not (yet) have a slash-command trigger for an on-demand full run.

Skips the handful of tests that dominate runtime on a genome-scale model
(flux-variability / loopless-MILP consistency checks, and the matrix-rank /
null-space tests, which are O(n^3)) so a pull request gets a score within a
normal CI budget.

Writes the total score and section scores to ``<out>/memote_score.md``
(diff-friendly markdown, read back by ``build_qc_report.py``) and the full
scored result to ``memote_result.json`` next to it -- not committed (see
model-qc.yml, which uploads it as a build artifact instead) to avoid
bloating the repository with a JSON blob on every pull request.

Usage:
    python code/python/tests/ci/memote_snapshot.py --out data/testResults
"""
from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

import cobra
import memote.suite.api as api
from memote.suite.reporting import ReportConfiguration, SnapshotReport
from raven_toolbox.io import read_yaml_model

from yeastgem import load_yeast_yaml

# The tests that dominate MEMOTE runtime on a genome-scale model:
#  * consistency: MILP / flux-variability / per-metabolite optimisation over
#    the whole model (stoichiometric consistency, energy cycles, blocked
#    reactions, open-bound producibility, ...).
#  * matrix: rank / null-space of the stoichiometric matrix, O(n^3) and
#    intractable at genome scale.
# Skipping these keeps a pull-request run within a normal CI budget.
_SLOW_TESTS = [
    # consistency (test_consistency.py)
    "test_stoichiometric_consistency",
    "test_unconserved_metabolites",
    "test_inconsistent_min_stoichiometry",
    "test_detect_energy_generating_cycles",
    "test_find_stoichiometrically_balanced_cycles",
    "test_blocked_reactions",
    "test_find_reactions_unbounded_flux_default_condition",
    "test_find_metabolites_not_produced_with_open_bounds",
    "test_find_metabolites_not_consumed_with_open_bounds",
    # matrix (test_matrix.py)
    "test_number_independent_conservation_relations",
    "test_matrix_rank",
    "test_degrees_of_freedom",
]


def _total_score(scored: dict) -> float | None:
    """Pull the total score (0-1) out of a scored MEMOTE result."""
    score = scored.get("score")
    if isinstance(score, dict) and isinstance(score.get("total_score"), (int, float)):
        return float(score["total_score"])
    if isinstance(scored.get("total_score"), (int, float)):
        return float(scored["total_score"])
    return None


def _section_rows(scored: dict) -> list[tuple[str, float]]:
    """Per-section scores for the summary table, best-effort."""
    score = scored.get("score")
    sections = score.get("sections") if isinstance(score, dict) else None
    rows = []
    for section in sections or []:
        name = section.get("section") or section.get("title")
        value = section.get("score")
        if name is not None and isinstance(value, (int, float)):
            rows.append((str(name), float(value)))
    return rows


def _test_metric(scored: dict, test_id: str) -> float | None:
    """The 0-1 metric of a single MEMOTE test (parametrised tests averaged)."""
    test = (scored.get("tests") or {}).get(test_id)
    if not isinstance(test, dict):
        return None
    metric = test.get("metric")
    if isinstance(metric, (int, float)):
        return float(metric)
    if isinstance(metric, dict):
        values = [v for v in metric.values() if isinstance(v, (int, float))]
        return sum(values) / len(values) if values else None
    return None


def _test_title(scored: dict, test_id: str) -> str:
    test = (scored.get("tests") or {}).get(test_id) or {}
    return str(test.get("title") or test_id)


def _detailed_rows(scored: dict, config: ReportConfiguration
                    ) -> list[tuple[str, str, float]]:
    """Per-test scores grouped by section: ``(section, test, metric)``.

    Uses the scoring configuration's section -> cases mapping to place
    each test, and the per-test metric from the scored result. Falls back
    to a flat list of every scored test if that mapping is not where
    expected, so the detail is still available even if less tidy.
    """
    try:
        sections = config["cards"]["scored"]["sections"]
    except (KeyError, TypeError, AttributeError):
        sections = None
    rows: list[tuple[str, str, float]] = []
    if isinstance(sections, dict):
        for section_id, section in sections.items():
            if not isinstance(section, dict):
                continue
            title = str(section.get("title") or section_id)
            for case in section.get("cases") or []:
                metric = _test_metric(scored, case)
                if metric is not None:
                    rows.append((title, _test_title(scored, case), metric))
    if rows:
        return rows
    tests = scored.get("tests")
    if isinstance(tests, dict):
        for test_id in tests:
            metric = _test_metric(scored, test_id)
            if metric is not None:
                rows.append(("", _test_title(scored, test_id), metric))
    return rows


def snapshot(model: cobra.Model) -> tuple[dict, ReportConfiguration]:
    """Run the MEMOTE fast subset on ``model`` and return the scored result."""
    with tempfile.TemporaryDirectory() as tmp:
        sbml_path = str(Path(tmp) / "model.xml")
        cobra.io.write_sbml_model(model, sbml_path)
        model_obj, sbml_ver, _ = api.validate_model(sbml_path)
        print(f"Running MEMOTE (fast subset, skipping {len(_SLOW_TESTS)} slow "
              "tests)", flush=True)
        _, result = api.test_model(
            model_obj, sbml_version=sbml_ver, results=True,
            skip=_SLOW_TESTS, solver_timeout=120,
        )
    config = ReportConfiguration.load()
    scored = SnapshotReport(result=result, configuration=config).result
    return scored, config


def render(scored: dict, config: ReportConfiguration) -> str:
    total = _total_score(scored)
    if total is None:
        return "# MEMOTE snapshot\n\nTotal score: unavailable (see workflow log).\n"
    pct = total * 100 if total <= 1 else total
    lines = [
        "# MEMOTE snapshot",
        "",
        f"Skipped (slow) tests: {', '.join(_SLOW_TESTS)}.",
        "",
        f"**Total score: {pct:.1f}%**",
    ]
    rows = _section_rows(scored)
    if rows:
        lines += ["", "| Section | Score |", "| --- | ---: |"]
        lines += [f"| {name} | {value * 100:.1f}% |" for name, value in rows]

    # A plain table here, not the <details> HTML the pull-request comment
    # wraps it in: this file is the source of data, build_qc_report.py
    # decides how to present it (and re-parses this same table with the
    # same 3-column pattern for the collapsible block).
    try:
        detailed = _detailed_rows(scored, config)
    except Exception as exc:  # detail is optional, never fail the run on it
        print(f"::warning::Could not build the detailed MEMOTE scores ({exc}).")
        detailed = []
    if detailed:
        lines += ["", "### Detailed scores", "", "| Section | Test | Score |",
                  "| --- | --- | ---: |"]
        lines += [f"| {section} | {test} | {metric * 100:.1f}% |"
                  for section, test, metric in detailed]
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=None,
                        help="model to score (default: this checkout's)")
    parser.add_argument("--out", type=Path, default=Path("data/testResults"),
                        help="directory to write memote_score.md into")
    args = parser.parse_args()

    model = (
        cobra.io.read_sbml_model(str(args.model))
        if args.model is not None and args.model.suffix == ".xml"
        # raven_toolbox.io.read_yaml_model, not cobra.io.load_yaml_model:
        # the latter silently drops RAVEN-only fields (deltaG,
        # confidence_score, notes, inchis, model.id/name/version).
        else read_yaml_model(str(args.model))
        if args.model is not None
        else load_yeast_yaml()
    )

    scored, config = snapshot(model)

    args.out.mkdir(parents=True, exist_ok=True)
    with open(args.out / "memote_result.json", "w", encoding="utf-8") as fh:
        json.dump(scored, fh, default=str)

    summary = render(scored, config)
    (args.out / "memote_score.md").write_text(summary, encoding="utf-8")
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
