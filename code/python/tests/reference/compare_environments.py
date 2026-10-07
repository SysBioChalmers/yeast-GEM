"""Compare MATLAB applyEnvironment with Python yeastgem.conditions.apply.

Usage: python compare_environments.py <dir>
where <dir> holds the output of dumpEnvironments.m. Every environment in
data/conditions must have a dump, and bounds and stoichiometry must be
identical (coefficients within 1e-9). Exit code 0 if all are identical.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

from yeastgem import conditions, load_yeast_yaml


def main(dump_dir: Path) -> int:
    names = sorted(p.stem for p in (Path(conditions.__file__).resolve().parents[3]
                                    / "data" / "conditions").glob("*.yml"))
    base = load_yeast_yaml()
    all_same = True
    for name in names:
        bounds_file = dump_dir / f"{name}_bounds.tsv"
        if not bounds_file.is_file():
            print(f"{name:22s} no MATLAB dump")
            all_same = False
            continue
        model = conditions.apply(base.copy(), name)
        with open(bounds_file, newline="") as handle:
            ml_bounds = {r["rxn"]: (float(r["lb"]), float(r["ub"]))
                         for r in csv.DictReader(handle, delimiter="\t")}
        with open(dump_dir / f"{name}_S.tsv", newline="") as handle:
            ml_s = {(r["met"], r["rxn"]): float(r["coef"])
                    for r in csv.DictReader(handle, delimiter="\t")}
        py_bounds = {r.id: r.bounds for r in model.reactions}
        py_s = {(m.id, r.id): c for r in model.reactions for m, c in r.metabolites.items()}
        bound_diff = [k for k in ml_bounds.keys() | py_bounds.keys()
                      if ml_bounds.get(k) != py_bounds.get(k)]
        s_diff = [k for k in ml_s.keys() | py_s.keys()
                  if abs(ml_s.get(k, 0.0) - py_s.get(k, 0.0)) > 1e-9]
        same = not bound_diff and not s_diff
        all_same &= same
        print(f"{name:22s} bounds differ: {len(bound_diff):4d}  S differs: {len(s_diff):4d}"
              f"  {'identical' if same else 'DIFFERENT'}")
    return 0 if all_same else 1


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1])))
