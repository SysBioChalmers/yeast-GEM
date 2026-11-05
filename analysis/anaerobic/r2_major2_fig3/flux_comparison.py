"""Statistics of the 13C flux comparison (Figure 3C-D).

Reads repro/fig3cd_points.tsv (measured and simulated flux, % of the glucose
uptake, for every measured flux of each data set) and writes
results/flux_comparison.tsv: for each model the number of fluxes, the mean
relative error of the large fluxes (measured above 50) and of the small fluxes
(measured below 5), and the largest absolute deviations. Glucose uptake is
imposed and excluded.

Usage: python flux_comparison.py [results_dir]
"""
import csv
import os
import sys
from pathlib import Path

RESULTS = Path(sys.argv[1] if len(sys.argv) > 1 else os.environ.get("ANAEROBIC_RESULTS", Path(__file__).resolve().parents[1] / "results"))
GLUCOSE = "r_1714"

rows = list(csv.DictReader(open(RESULTS / "repro" / "fig3cd_points.tsv", encoding="utf-8"), delimiter="\t"))
out = ["variant\tn_fluxes\tn_large\tMRE_large\tn_small\tMRE_small\tlargest_absolute_deviations"]
for variant in dict.fromkeys(r["variant"] for r in rows):
    pts = [(r["dataset"], r["rxn"], float(r["measured"]), float(r["simulated"]))
           for r in rows if r["variant"] == variant and r["rxn"] != GLUCOSE]
    mre = lambda sel: sum(abs(s - m) / abs(m) for _, _, m, s in sel) / len(sel) if sel else float("nan")
    large = [p for p in pts if abs(p[2]) > 50]
    small = [p for p in pts if 0 < abs(p[2]) < 5]
    worst = sorted(pts, key=lambda p: -abs(p[3] - p[2]))[:3]
    text = "; ".join(f"{d} {r} {s - m:+.1f}" for d, r, m, s in worst)
    out.append(f"{variant}\t{len(pts)}\t{len(large)}\t{mre(large):.3f}\t{len(small)}\t{mre(small):.3f}\t{text}")
(RESULTS / "flux_comparison.tsv").write_text("\n".join(out) + "\n", encoding="utf-8")
print("\n".join(out))
