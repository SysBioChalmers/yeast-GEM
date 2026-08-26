"""Level-2 parity gate — Python validation metrics match the committed
MATLAB-produced results.

Computes growth R², essential-gene accuracy and confusion matrix, and
anaerobic flux R², then checks each against the values published in
``data/testResults/README.md`` within tolerance. That file is the
human-readable record of the model tests, regenerated at release time,
so there is a single set of committed numbers rather than a separate
machine-readable copy that can drift away from it.

Run locally:
    python code/python/tests/ci/check_metrics.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from yeastgem import conditions, model_tests, read_yeast_model
from yeastgem.io import REPO_PATH

_DEFAULT_RESULTS = REPO_PATH / "data" / "testResults" / "README.md"

# Tolerances live here rather than beside the results, because they are a
# property of the comparison and not of the model.
#
# Two independent sources of drift are absorbed:
#
# 1. R² metrics differ at ~1e-4 between solvers. The reference is produced
#    by MATLAB with Gurobi; CI resolves Gurobi when the organisation licence
#    is available and GLPK otherwise.
# 2. The essential-gene confusion matrix differs from the MATLAB reference by
#    one gene (Python 933/65/94/15 vs MATLAB 934/65/94/14). This one is *not*
#    solver drift: it is unchanged when Python also uses Gurobi, so it
#    reflects a difference between cobrapy's single_gene_deletion and RAVEN's
#    findGeneDeletions around the 1e-6 growth-ratio threshold. Tracked
#    separately; gene_count_abs absorbs it in the meantime.
_TOL_R2 = 5e-3
_TOL_ACCURACY = 5e-3
_TOL_GENE_COUNT = 2

# Row labels in the results table, mapped to the keys used below.
_ROWS = {
    "growth prediction r2": "growth_r2",
    "anaerobic flux prediction r2": "anaerobic_flux_r2",
    "gene essentiality accuracy": "accuracy",
    "true non-essential genes": "tp",
    "true essential genes": "tn",
    "false non-essential genes": "fp",
    "false essential genes": "fn",
}


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


def main(results_path: Path = _DEFAULT_RESULTS) -> int:
    ref = parse_results(results_path)

    print(f"Reference: {results_path}")
    model = read_yeast_model()
    # No solver is pinned anywhere: cobrapy takes whatever optlang resolves,
    # which is Gurobi when the licence is available and GLPK otherwise. Since
    # that choice is what the tolerances above are absorbing, report it rather
    # than leaving it to be inferred.
    print(f"Solver: {model.solver.interface.__name__}")

    print("Computing growth R² ...")
    growth_r2 = model_tests.growth(model.copy())
    print(f"  Python: {growth_r2:.6g}  Reference: {ref['growth_r2']:.6g}")

    print("Computing essential_genes ...")
    result = model_tests.essential_genes(model.copy())
    print(f"  Python accuracy: {result.accuracy:.6g}  Reference: {ref['accuracy']:.6g}")
    print(f"  Python TP/TN/FP/FN: "
          f"{len(result.tp)}/{len(result.tn)}/{len(result.fp)}/{len(result.fn)}  "
          f"Reference: {ref['tp']:.0f}/{ref['tn']:.0f}/{ref['fp']:.0f}/{ref['fn']:.0f}")

    print("Computing anaerobic_flux_predictions ...")
    anaerobic = model.copy()
    conditions.apply(anaerobic, "anaerobic")
    af_r2, _af_mre = model_tests.anaerobic_flux_predictions(anaerobic)
    print(f"  Python: {af_r2:.6g}  Reference: {ref['anaerobic_flux_r2']:.6g}")

    checks: list[tuple[str, float, float, float]] = [
        ("growth R2", growth_r2, ref["growth_r2"], _TOL_R2),
        ("anaerobic flux R2", af_r2, ref["anaerobic_flux_r2"], _TOL_R2),
        ("gene essentiality accuracy", result.accuracy, ref["accuracy"], _TOL_ACCURACY),
        ("true non-essential genes", len(result.tp), ref["tp"], _TOL_GENE_COUNT),
        ("true essential genes", len(result.tn), ref["tn"], _TOL_GENE_COUNT),
        ("false non-essential genes", len(result.fp), ref["fp"], _TOL_GENE_COUNT),
        ("false essential genes", len(result.fn), ref["fn"], _TOL_GENE_COUNT),
    ]

    failures: list[str] = []
    for name, actual, expected, tol_abs in checks:
        diff = abs(actual - expected)
        if diff > tol_abs:
            failures.append(
                f"  FAIL {name}: |{actual:.6g} - {expected:.6g}| = "
                f"{diff:.6g}  > tol {tol_abs:.6g}"
            )

    if failures:
        print("\nMetric parity FAILED:")
        for msg in failures:
            print(msg)
        return 1
    print("\nAll metric-parity checks passed.")
    return 0


if __name__ == "__main__":
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else _DEFAULT_RESULTS
    raise SystemExit(main(path))
