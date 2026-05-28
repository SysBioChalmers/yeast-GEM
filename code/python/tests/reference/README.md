# Reference bundle — MATLAB-produced fixtures

This directory holds the **MATLAB-produced reference artifacts** that the
Python toolchain is compared against by the CI level-1 (semantic
equality) and level-2 (metric parity) gates.

## Contents (when populated)

- `yeast-GEM.xml` — the committed model saved by the MATLAB toolchain
  (RAVEN + `commitYeastModel`). Identical content, MATLAB-authored
  formatting; the comparator ignores formatting differences but checks
  semantic equality byte-aware.
- `metrics.yml` — reference values for the level-2 gate:
  ```yaml
  aerobic_growth: 0.0876...        # objective at optimum, aerobic minimal
  chemostat_r2: 0.97...            # growth.m R² across 4 conditions
  essential_genes:
    tp: 119
    tn: 922
    fp: 31
    fn: 65
  anaerobic_flux_r2: 0.95...
  biomass_fractions:
    X: 1.0
    P: 0.46
    C: 0.31
    R: 0.06
    D: 0.005
    L: 0.095
    I: 0.025
    F: 0.0045
  ```
- `provenance.yml` — MATLAB / RAVEN / solver versions and the git SHA
  the artifacts were generated from.

## Regeneration

The reference bundle is **not produced per-PR**. It is regenerated at
the start of each release cycle (or whenever a behavior change in the
committed model forces a refresh) by running
[`regenerate.m`](regenerate.m) in MATLAB with RAVEN on the same git
SHA the model is committed at.

```matlab
cd code/python/tests/reference
regenerate
```

The resulting `yeast-GEM.xml`, `metrics.yml` and `provenance.yml` are
committed to this directory in the same PR that introduced the change.

## CI usage

- The `matlab-reference-compare` job in
  [.github/workflows/python.yml](../../../../.github/workflows/python.yml)
  is currently `if: false` (gate disabled until the bundle is first
  seeded). Once `yeast-GEM.xml` is committed here, flip that gate to
  `true` and the level-1 comparison becomes a required check.
- The level-2 metric-parity gate (not yet wired) will load
  `metrics.yml` and compare Python-computed metrics to the reference
  values within the tolerances defined in
  [PORTING_PLAN.md](../../PORTING_PLAN.md).

## Why this is in MATLAB, not Python

Per the lock-step parity policy, the MATLAB toolchain is the
canonical source for the committed artifact during the transition. The
reference bundle locks in "this is what MATLAB produces" so the Python
port can be validated independently. When both toolchains have full
parity and a single-language production owner is chosen, this
direction may flip — until then, MATLAB seeds, Python verifies.
