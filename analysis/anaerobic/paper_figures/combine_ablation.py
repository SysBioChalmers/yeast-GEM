"""Combine the files <prefix>_<n>.tsv (parallel jobs of curationAblation) into
<prefix>.tsv. Usage: python combine_ablation.py [results_dir] [prefix]
prefix: curation_ablation (default) or curation_pairs."""
import csv
import os
import sys
from pathlib import Path

R = Path(sys.argv[1] if len(sys.argv) > 1 else os.environ.get("ANAEROBIC_RESULTS", Path(__file__).resolve().parents[1] / "results"))
PREFIX = sys.argv[2] if len(sys.argv) > 2 else "curation_ablation"
files = sorted(R.glob(PREFIX + "_*.tsv"), key=lambda p: int(p.stem.rsplit("_", 1)[1]))
rows = [r for f in files for r in csv.DictReader(open(f, encoding="utf-8"), delimiter="\t")]
with open(R / (PREFIX + ".tsv"), "w", encoding="utf-8", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
print(len(rows), "models")
