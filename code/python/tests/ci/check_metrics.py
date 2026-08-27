"""Validation metrics, compared against the target branch.

Computes growth R², gene-essentiality accuracy and its confusion matrix,
anaerobic flux R² and the anaerobic exchange metrics, for this branch's
model and for the target branch's model, and reports the difference.

Comparing against the target branch rather than against committed
reference numbers matters here. ``data/testResults/README.md`` records
the values as of the last *release*, so on any curation pull request the
model has legitimately moved away from it and the comparison degenerates
into "is the tolerance wide enough". Measuring both sides in the same run
-- same code, same reference data, same solver -- isolates what this
branch actually changed.

The reference data (chemostat measurements, the deletion collection, the
anaerobic exchange rates) always comes from the current checkout: only
the model differs between the two columns.

With no ``--base-model`` the metrics are checked against the committed
table instead, which is what the release routine wants.

Run locally:
    python code/python/tests/ci/check_metrics.py
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import cobra
import matplotlib

matplotlib.use("Agg")

from yeastgem import conditions, model_tests, read_yeast_model
from yeastgem.io import REPO_PATH

_DEFAULT_RESULTS = REPO_PATH / "data" / "testResults" / "README.md"

# Tolerances are a property of the comparison, not of the model, so they
# live here rather than beside the results.
#
# Against the target branch they absorb solver noise only, since both sides
# run in the same job. Against the committed table they additionally absorb
# the reference having been produced by MATLAB with Gurobi, and a one-gene
# difference in the confusion matrix that is not solver drift: it persists
# when Python also uses Gurobi, and reflects cobrapy's single_gene_deletion
# differing from RAVEN's findGeneDeletions around the 1e-6 growth-ratio
# threshold.
_TOL_R2 = 5e-3
_TOL_ACCURACY = 5e-3
_TOL_GENE_COUNT = 2
_TOL_RELATIVE_ERROR = 1e-2
_TOL_RATIO = 5e-2

# (key, comment label, direction, tolerance)
#
# "direction" says which way is an improvement, and is what makes a delta
# readable: a rising R² is good, a rising error is not, and a rising
# false-negative count is not. A move within tolerance is reported as
# unchanged.
_METRICS = [
    ("growth_r2", "Growth prediction R2", "higher", _TOL_R2),
    ("anaerobic_flux_r2", "Anaerobic flux prediction R2", "higher", _TOL_R2),
    ("exchange_mean_relative_error", "Anaerobic exchange mean relative error",
     "lower", _TOL_RELATIVE_ERROR),
    ("exchange_within_error", "Anaerobic exchange within error", "higher",
     _TOL_RELATIVE_ERROR),
    ("ammonium_per_atpase", "Ammonium per ATPase", "none", _TOL_RATIO),
    ("accuracy", "Gene essentiality accuracy", "higher", _TOL_ACCURACY),
    ("tp", "True non-essential genes", "higher", _TOL_GENE_COUNT),
    ("tn", "True essential genes", "higher", _TOL_GENE_COUNT),
    ("fp", "False non-essential genes", "lower", _TOL_GENE_COUNT),
    ("fn", "False essential genes", "lower", _TOL_GENE_COUNT),
]

# Row labels in the committed table, mapped to the keys used here.
_ROWS = {label.lower(): key for key, label, _dir, _tol in _METRICS}


def parse_results(path: Path) -> dict[str, float]:
    """Read the metric table out of data/testResults/README.md.

    Expects rows of the form ``| Metric name | value |``. Anything else in
    the file is ignored, so the prose around the table is free to change.
    """
    found: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        match = re.match(r"^\s*\|\s*(.+?)\s*\|\s*([-\d.eE+]+)\s*\|\s*$", line)
        if not match:
            continue
        key = _ROWS.get(match.group(1).strip().lower())
        if key is not None:
            found[key] = float(match.group(2))

    missing = sorted(set(_ROWS.values()) - set(found))
    if missing:
        raise SystemExit(
            f"{path} is missing rows for: {', '.join(missing)}. "
            "The table must keep one '| name | value |' row per metric."
        )
    return found


def compute(model) -> dict[str, float]:
    """Every validation metric for one model."""
    growth_r2 = model_tests.growth(model.copy())
    essential = model_tests.essential_genes(model.copy())

    anaerobic = model.copy()
    conditions.apply(anaerobic, "anaerobic")
    anaerobic_r2, _mre = model_tests.anaerobic_flux_predictions(anaerobic)
    exchange = model_tests.plot_anaerobic(anaerobic, plot=False)

    return {
        "growth_r2": float(growth_r2),
        "anaerobic_flux_r2": float(anaerobic_r2),
        "exchange_mean_relative_error": exchange.mean_relative_error,
        "exchange_within_error": exchange.fraction_within_error,
        "ammonium_per_atpase": exchange.ammonium_per_atpase,
        "accuracy": float(essential.accuracy),
        "tp": float(len(essential.tp)),
        "tn": float(len(essential.tn)),
        "fp": float(len(essential.fp)),
        "fn": float(len(essential.fn)),
    }


def _verdict(value: float, base: float, direction: str, tol: float):
    """Return ``(delta_text, icon, is_regression)`` for one metric."""
    change = value - base
    if abs(change) <= tol:
        return "0", ":white_check_mark:", False
    text = f"{change:+.4g}"
    if direction == "none":
        return text, ":warning:", False
    improved = (change > 0) if direction == "higher" else (change < 0)
    return text, (":white_check_mark:" if improved else ":x:"), not improved


def render(current: dict, base: dict, base_ref: str, solver: str,
           run_url: str = "") -> str:
    """Render the comparison for a pull-request comment."""
    rows, regressions, moved = [], 0, 0
    for key, label, direction, tol in _METRICS:
        value = current[key]
        if key not in base:
            rows.append(f"| {label} | {value:.6g} | — | new | :grey_question: |")
            continue
        delta, icon, regression = _verdict(value, base[key], direction, tol)
        regressions += regression
        moved += delta != "0"
        rows.append(
            f"| {label} | {value:.6g} | {base[key]:.6g} | {delta} | {icon} |"
        )

    if regressions:
        verdict = (
            f":x: **{regressions} metric(s) got worse vs `{base_ref}`.** "
            "Review the :x: rows."
        )
    elif moved:
        verdict = (
            f":white_check_mark: **No regressions vs `{base_ref}`** "
            f"({moved} metric(s) moved, none for the worse)."
        )
    else:
        verdict = (
            f":white_check_mark: **Unchanged vs `{base_ref}`** — every metric "
            "is within tolerance."
        )

    lines = [
        "## Validation metrics",
        "",
        verdict,
        "",
        f"| Metric | This branch | `{base_ref}` | &Delta; | |",
        "| --- | ---: | ---: | ---: | :---: |",
        *rows,
        "",
        f"Solver: `{solver}`. A change within tolerance is reported as 0. "
        "Direction matters: a rising R2 or accuracy is an improvement, a "
        "rising error or false-call count is not.",
        "",
        "_Both columns are computed in this run, against the same reference "
        "data — only the model differs._",
    ]
    if run_url:
        lines += ["", f"[Full workflow run]({run_url})"]
    return "\n".join(lines) + "\n"


def _load(path: Path | None):
    if path is None:
        return read_yeast_model()
    if path.suffix in {".yml", ".yaml"}:
        return cobra.io.load_yaml_model(str(path))
    return cobra.io.read_sbml_model(str(path))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=None,
                        help="model to measure (default: this checkout's)")
    parser.add_argument("--base-model", type=Path, default=None,
                        help="target-branch model; when given, the comparison "
                             "is against this instead of the committed table")
    parser.add_argument("--results", type=Path, default=_DEFAULT_RESULTS,
                        help="committed results table, used when there is no "
                             "--base-model")
    parser.add_argument("--emit", type=Path, default=None,
                        help="write the measured values as key<TAB>value, for "
                             "build_qc_report.py to compare and render")
    parser.add_argument("--markdown", type=Path, default=None,
                        help="write the comparison here, for a PR comment")
    parser.add_argument("--base-ref", default="the target branch")
    parser.add_argument("--run-url", default="")
    args = parser.parse_args()

    model = _load(args.model)
    solver = model.solver.interface.__name__
    print(f"Solver: {solver}")

    print("Measuring this branch ...")
    current = compute(model)

    if args.emit is not None:
        args.emit.parent.mkdir(parents=True, exist_ok=True)
        args.emit.write_text(
            "".join(f"{k}\t{current[k]}\n" for k in sorted(current)),
            encoding="utf-8",
        )
        print(f"Wrote {args.emit}")
        # Emitting is for the report, which does its own comparison against
        # the target branch. Comparing here as well would mean measuring the
        # base model twice.
        return 0

    if args.base_model is not None:
        print(f"Measuring {args.base_ref} ...")
        base = compute(_load(args.base_model))
        source = args.base_ref
    else:
        base = parse_results(args.results)
        source = str(args.results)

    width = max(len(label) for _k, label, _d, _t in _METRICS)
    regressions = []
    for key, label, direction, tol in _METRICS:
        if key not in base:
            continue
        delta, _icon, regression = _verdict(current[key], base[key], direction, tol)
        print(f"  {label:<{width}}  {current[key]:>12.6g}  vs "
              f"{base[key]:>12.6g}   {delta}")
        if regression:
            regressions.append(label)

    if args.markdown is not None:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(
            render(current, base, source, solver, args.run_url),
            encoding="utf-8",
        )
        print(f"\nWrote summary to {args.markdown}")

    if regressions:
        print(f"\nFAILED — worse than {source}: {', '.join(regressions)}")
        return 1
    print(f"\nNo regressions against {source}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
