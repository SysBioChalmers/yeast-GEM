"""Convergence diagnostics over independent RHMC chains (seeds 1-4).

For every (variant, scenario) with chain files in results/sampling (seed 1 =
<variant>__<scenario>__opt<f>.mat, seeds 2-4 = ...__seed<n>.mat), on the
cycle-free samples W:
- rank-normalised split R-hat and bulk/tail ESS (Vehtari et al 2021, ArviZ),
  for all reactions that vary and for the reactions the manuscript reports;
- per-chain 2.5/50/97.5% quantiles of reported fluxes;
- for 13C scenarios, mean CRPS per chain and pooled (fluxes as |% of glucose
  uptake|, as analyzeSampling.m), with and without the glucose points.
Writes results/chains/{diagnostics,quantiles,crps}.tsv.

Usage: python chain_diagnostics.py [optFrac]
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import arviz as az
import h5py
import numpy as np

STUDY = Path(__file__).resolve().parents[1]
SAMPLES = STUDY / "results" / "sampling"
OUT = STUDY / "results" / "chains"
FLUX = STUDY.parents[1] / "data" / "physiology" / "flux_data_anaerobic.tsv"
REPORTED = {"r_4041": "growth", "r_1808": "glycerol", "r_1761": "ethanol", "r_1672": "CO2",
            "r_2056": "succinate", "r_0227": "PM H+-ATPase", "r_2131": "Idp1", "r_0658": "Idh1/2",
            "r_1810": "glycine export", "r_0226": "mito ATP synthase", "r_1832": "H+ export"}


def load(path: Path) -> tuple[np.ndarray, list[str]]:
    with h5py.File(path, "r") as f:
        W = np.asarray(f["W"])                      # draws x rxns (MATLAB rxns x draws)
        refs = f["rxns"][0]
        rxns = ["".join(chr(c) for c in np.asarray(f[r]).ravel()) for r in refs]
    return W, rxns


def crps(x: np.ndarray, y: float) -> float:
    x = np.sort(x)
    n = len(x)
    return float(np.mean(np.abs(x - y)) - np.sum((2 * np.arange(1, n + 1) - n - 1) * x) / n**2)


def flux_data() -> dict[str, list[tuple[str, float]]]:
    data: dict[str, list[tuple[str, float]]] = {}
    with open(FLUX, newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) > 6 and row[6].startswith("r_"):
                data.setdefault(row[5], []).append((row[6], abs(float(row[4]))))
    return data


def main(opt: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    meas = flux_data()
    diag_rows, q_rows, c_rows = [], [], []
    bases = sorted(p for p in SAMPLES.glob(f"*__opt{opt}.mat"))
    for base in bases:
        variant, scenario = base.name.split("__")[:2]
        files = [base] + [SAMPLES / f"{variant}__{scenario}__opt{opt}__seed{s}.mat" for s in (2, 3, 4)]
        files = [f for f in files if f.exists()]
        if len(files) < 2:
            continue
        chains, rxns = [], None
        for f in files:
            try:
                W, r = load(f)
            except (KeyError, OSError):
                break                               # infeasible scenario: no samples
            if rxns is not None and r != rxns:
                raise ValueError(f"reaction order differs in {f}")
            rxns = r
            chains.append(W)
        if len(chains) < 2:
            continue
        n = min(c.shape[0] for c in chains)
        arr = np.stack([c[:n] for c in chains])     # chain x draw x rxn
        span = arr.max(axis=(0, 1)) - arr.min(axis=(0, 1))
        moving = span > 1e-9
        idx = {r: i for i, r in enumerate(rxns)}
        reported = [idx[r] for r in REPORTED if r in idx]
        compared = [idx[r] for r, _ in meas.get(scenario.removeprefix("flux_"), []) if r in idx]
        ds = az.convert_to_dataset({"v": arr[:, :, moving]})
        rhat = az.rhat(ds, method="rank")["v"].values
        bulk = az.ess(ds, method="bulk")["v"].values
        tail = az.ess(ds, method="tail")["v"].values
        mov_idx = np.flatnonzero(moving)
        pos = {j: k for k, j in enumerate(mov_idx)}
        for label, sel in (("all varying", list(mov_idx)),
                           ("reported", [j for j in reported if j in pos]),
                           ("compared 13C", [j for j in compared if j in pos])):
            if not sel:
                continue
            k = [pos[j] for j in sel]
            diag_rows.append({"variant": variant, "scenario": scenario, "subset": label, "chains": len(chains),
                              "draws_per_chain": n, "n_rxns": len(k),
                              "rhat_max": round(float(np.nanmax(rhat[k])), 4),
                              "frac_rhat_gt_1.01": round(float(np.mean(rhat[k] > 1.01)), 4),
                              "ess_bulk_min": round(float(np.nanmin(bulk[k]))),
                              "ess_tail_min": round(float(np.nanmin(tail[k])))})
        for r, name in REPORTED.items():
            if r not in idx:
                continue
            j = idx[r]
            for c in range(len(chains)):
                q = np.percentile(arr[c, :, j], [2.5, 50, 97.5])
                q_rows.append({"variant": variant, "scenario": scenario, "rxn": r, "name": name, "chain": c + 1,
                               "q2_5": q[0], "median": q[1], "q97_5": q[2]})
        if scenario.startswith("flux_"):
            ds_name = scenario.removeprefix("flux_")
            g = idx["r_1714"]
            for c in list(range(len(chains))) + ["pooled"]:
                block = arr[c] if c != "pooled" else arr.reshape(-1, arr.shape[2])
                vals = {}
                for r, y in meas.get(ds_name, []):
                    if r in idx:
                        x = np.abs(100 * block[:, idx[r]] / block[:, g])
                        vals[r] = crps(x, y)
                no_glc = [v for r, v in vals.items() if r != "r_1714"]
                c_rows.append({"variant": variant, "dataset": ds_name, "chain": c if c == "pooled" else c + 1,
                               "n": len(vals), "crps_mean": np.mean(list(vals.values())),
                               "crps_mean_no_glucose": np.mean(no_glc)})
        print(f"{variant:10s} {scenario:20s} chains {len(chains)} draws {n}: "
              + "; ".join(f"{d['subset']} R-hat max {d['rhat_max']} ESS bulk/tail min {d['ess_bulk_min']}/{d['ess_tail_min']}"
                          for d in diag_rows if d["variant"] == variant and d["scenario"] == scenario))
    for name, rows in (("diagnostics", diag_rows), ("quantiles", q_rows), ("crps", c_rows)):
        if rows:
            with open(OUT / f"{name}_opt{opt}.tsv", "w", newline="") as handle:
                w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
                w.writeheader()
                w.writerows(rows)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "0.99999")
