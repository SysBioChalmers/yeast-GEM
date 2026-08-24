"""Level-2 parity gate — Python validation metrics match the
committed MATLAB-produced reference.

Computes growth R², essential-gene accuracy / sensitivity /
specificity / MCC + confusion matrix, and anaerobic flux R², then
checks each value against the committed reference (see
``code/python/tests/reference/metrics.json``) within tolerance.

Run locally:
    python code/python/tests/ci/check_metrics.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from yeastgem import conditions, model_tests, read_yeast_model

_DEFAULT_REFERENCE = Path(__file__).resolve().parents[1] / "reference" / "metrics.json"


def main(reference_path: Path = _DEFAULT_REFERENCE) -> int:
    ref = json.loads(reference_path.read_text())
    tol = ref["tolerances"]

    print(f"Reference: {reference_path} (source={ref.get('source_commit', '?')})")
    model = read_yeast_model()
    # No solver is pinned anywhere: cobrapy takes whatever optlang resolves,
    # which is GLPK in CI and Gurobi on a machine with gurobipy installed.
    # Since that choice is what the tolerances below are absorbing, report it
    # rather than leaving it to be inferred.
    print(f"Solver: {model.solver.interface.__name__}")

    print("Computing growth R² ...")
    growth_r2 = model_tests.growth(model.copy())
    print(f"  Python: {growth_r2:.6g}  Reference: {ref['growth_r2']:.6g}")

    print("Computing essential_genes ...")
    result = model_tests.essential_genes(model.copy())
    eg_ref = ref["essential_genes"]
    print(f"  Python accuracy: {result.accuracy:.6g}  Reference: {eg_ref['accuracy']:.6g}")
    print(f"  Python TP/TN/FP/FN: "
          f"{len(result.tp)}/{len(result.tn)}/{len(result.fp)}/{len(result.fn)}  "
          f"Reference: {eg_ref['tp']}/{eg_ref['tn']}/{eg_ref['fp']}/{eg_ref['fn']}")

    print("Computing anaerobic_flux_predictions ...")
    anaerobic = model.copy()
    conditions.apply(anaerobic, "anaerobic")
    af_r2, _af_mre = model_tests.anaerobic_flux_predictions(anaerobic)
    print(f"  Python: {af_r2:.6g}  Reference: {ref['anaerobic_flux_r2']:.6g}")

    checks: list[tuple[str, float, float, float]] = [
        ("growth_r2", growth_r2, ref["growth_r2"], tol["r2_abs"]),
        ("essential_genes accuracy", result.accuracy, eg_ref["accuracy"], tol["accuracy_abs"]),
        ("essential_genes sensitivity", result.sensitivity,
         eg_ref["sensitivity_percent"], 1.0),
        ("essential_genes specificity", result.specificity,
         eg_ref["specificity_percent"], 1.0),
        ("essential_genes mcc", result.mcc, eg_ref["mcc"], tol["mcc_abs"]),
        ("essential_genes tp", len(result.tp), eg_ref["tp"], tol["gene_count_abs"]),
        ("essential_genes tn", len(result.tn), eg_ref["tn"], tol["gene_count_abs"]),
        ("essential_genes fp", len(result.fp), eg_ref["fp"], tol["gene_count_abs"]),
        ("essential_genes fn", len(result.fn), eg_ref["fn"], tol["gene_count_abs"]),
        ("anaerobic_flux_r2", af_r2, ref["anaerobic_flux_r2"], tol["r2_abs"]),
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
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else _DEFAULT_REFERENCE
    raise SystemExit(main(path))
