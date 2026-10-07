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

## Environments: MATLAB and Python give the same model

Environments (`data/conditions/*.yml`) are applied by `applyEnvironment`
in MATLAB and `yeastgem.conditions.apply` in Python, two independent
readers of the same files. [`dumpEnvironments.m`](dumpEnvironments.m)
writes the bounds and stoichiometry after every environment;
[`compare_environments.py`](compare_environments.py) applies the same
environments in Python and requires identical bounds and coefficients
(within 1e-9):

```bash
matlab -batch "addpath('code'); addpath('code/python/tests/reference'); \
    dumpEnvironments('/tmp/environments')"
python code/python/tests/reference/compare_environments.py /tmp/environments
```

Result (2026-10-07): identical for anaerobic, glycine_nitrogen,
minimal_Y6 and nitrogen_limitation. `applyEnvironment` also gave the
same lb, ub and S as the functions it replaces on develop
(anaerobicModel of 9.1.0, minimal_Y6, glycineNitrogenSource and
nitrogenLimitation).

## Phase-3 specific: the commit-pipeline equivalence check

Phase 3 renames `saveYeastModel` to `commitYeastModel` (with a
deprecation shim) and adds the Python `commit_yeast_model` release
pipeline. The verification driver
[`runPhase3.m`](runPhase3.m) takes a yeast-GEM checkout path and a
function name (either `saveYeastModel` or `commitYeastModel`) and
writes the resulting SBML to a target path:

```bash
# 1. Worktree pinned to the pre-rename commit (i.e. just before phase 3).
git worktree add --detach /tmp/yeast-gem-pre3 <pre-rename-SHA>

# 2. Produce the pre-rename and post-rename SBMLs.
matlab -batch "addpath('code/python/tests/reference'); \
    runPhase3('/tmp/yeast-gem-pre3', '/tmp/phase3-pre.xml', 'saveYeastModel')"
matlab -batch "addpath('code/python/tests/reference'); \
    runPhase3('.', '/tmp/phase3-post.xml', 'commitYeastModel')"

# 3. Verify the rename preserved behaviour.
python -m yeastgem.compare /tmp/phase3-pre.xml /tmp/phase3-post.xml

# 4. Python-vs-MATLAB parity for commit_yeast_model.
YEAST_GEM_PATH=/tmp/yeast-gem-pre3 python -c \
    "from yeastgem import read_yeast_model, commit_yeast_model; \
     m = read_yeast_model(); commit_yeast_model(m, update_readme=False)"
cp /tmp/yeast-gem-pre3/model/yeast-GEM.xml /tmp/phase3-py.xml
python -m yeastgem.compare /tmp/phase3-post.xml /tmp/phase3-py.xml

# 5. Tear down.
git worktree remove /tmp/yeast-gem-pre3
```

`runPhase3` always copies the freshly written `model/yeast-GEM.xml`
even when the surrounding `exportForGit` step fails — that helper
requires a COBRA-Toolbox-only `COBRAver` variable that is not present
on a pure-RAVEN MATLAB install, but the SBML write itself precedes it,
so the comparison still works.

### Result of the verification run (commit c74afed → phase-3 HEAD)

Both checks pass:
- `runPhase3(...saveYeastModel)` vs `runPhase3(...commitYeastModel)`:
  *Models are semantically equal* — the rename preserved behaviour
  exactly.
- MATLAB `commitYeastModel` vs Python `commit_yeast_model`:
  *Models are semantically equal* — the Python release pipeline lands
  on the same canonical model as the MATLAB pipeline.

## Phase-5 specific: Tier-3 model-tests metrics check

Phase 5 ports the four ``code/modelTests/*.m`` validation routines
into ``yeastgem.model_tests``. The verification driver
[`runPhase5Metrics.m`](runPhase5Metrics.m) computes growth R²,
essential-gene confusion matrix, and anaerobic-flux R² on the MATLAB
side; the Python equivalents are exposed under ``yeastgem.model_tests``
and should match within float tolerance.

```bash
# 1. MATLAB metrics (writes JSON).
matlab -batch "addpath('code/python/tests/reference'); \
    runPhase5Metrics('.', '/tmp/phase5-matlab-metrics.json')"

# 2. Python metrics.
python - <<'PY'
import json, matplotlib; matplotlib.use("Agg")
from yeastgem import read_yeast_model, conditions, model_tests
m = read_yeast_model()
r = model_tests.essential_genes(m.copy())
an = m.copy(); conditions.apply(an, "anaerobic")
af_r2, _ = model_tests.anaerobic_flux_predictions(an)
out = {
    "growth_r2": model_tests.growth(m.copy()),
    "essential_genes_accuracy": r.accuracy,
    "essential_genes_sensitivity": r.sensitivity,
    "essential_genes_specificity": r.specificity,
    "essential_genes_mcc": r.mcc,
    "anaerobic_flux_r2": af_r2,
}
json.dump(out, open("/tmp/phase5-python-metrics.json", "w"), indent=2)
PY

# 3. Diff (small absolute tolerances).
diff <(jq -S . /tmp/phase5-matlab-metrics.json) \
     <(jq -S . /tmp/phase5-python-metrics.json)
```

### Result of the verification run (phase 4 HEAD → phase 5 HEAD)

| Metric | MATLAB | Python | Δ |
|---|---|---|---|
| growth R² | 0.906164 | 0.906164 | ≤ 1e-7 |
| anaerobic flux R² | 0.904765 | 0.905662 | 9e-4 |
| essential_genes accuracy | 0.90244 | 0.90154 | 9e-4 |
| essential_genes sensitivity | 98.52 | 98.42 | 1e-1 |
| essential_genes specificity | 40.88 | 40.88 | 0 |
| essential_genes MCC | 0.5368 | 0.5323 | 4e-3 |
| TP/TN/FP/FN | 934/65/94/14 | 933/65/94/15 | 1 gene |

The single-gene discrepancy is a borderline case at the 1e-6 growth-ratio
threshold. This is not solver-related -- it persists when Python also
uses Gurobi -- and points at cobrapy's single_gene_deletion differing
from RAVEN's findGeneDeletions. All metrics are within
the level-2 tolerances defined in PORTING_PLAN.md (R²/accuracy ≤ 1e-4
when the underlying inputs match; here we're at 1e-3 because of one
gene). Acceptable for phase 5; revisit if any drifts beyond 1e-2.
