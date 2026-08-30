# yeastgem — Python interface to yeast-GEM

Python counterpart to the MATLAB code under [../](..). Builds on
[cobrapy](https://github.com/opencobra/cobrapy) and
[raven-toolbox](https://github.com/SysBioChalmers/raven-toolbox) — the
latter hosts the generic GEM utilities (model diffing, SBO term
assignment, condition / biomass / curation engines) that `yeastgem`
configures with the yeast-specific data files under
[../../data/](../../data/).

See [PORTING_PLAN.md](PORTING_PLAN.md) for the porting history and
[UPSTREAM_CANDIDATES.md](UPSTREAM_CANDIDATES.md) for the
function-level upstream tracking.

## Install (development)

```bash
pip install -e code/python/[dev]
```

The generic, organism-agnostic helpers live in
[raven-toolbox](https://github.com/SysBioChalmers/raven-python), which
is resolved from PyPI via the `raven-toolbox>=0.4.0` constraint in
[pyproject.toml](pyproject.toml).

## Quick start

```python
from yeastgem import read_yeast_model, commit_yeast_model

model = read_yeast_model()
print(model.optimize().objective_value)        # → ~0.088 / h on the default media

# Make some changes …
commit_yeast_model(model)                       # full release pipeline
```

`yeastgem` auto-locates the repo root via the package install path,
the `YEAST_GEM_PATH` environment variable, or a legacy `.env` file —
in that order. No additional setup needed for the common case.

## API map

| Area | Module | Highlights |
|---|---|---|
| **I/O** | [`yeastgem.io`](yeastgem/io.py) | `read_yeast_model`, `commit_yeast_model` (release pipeline: canonical state → SBML validity → aerobic + anaerobic growth → write `.xml` + ΔG CSVs). `write_yeast_model` is a deprecated forwarding shim. |
| **Comparison** | [`yeastgem.compare`](yeastgem/compare.py) | `compare_models` / `ComparisonReport` re-exported from `raven_toolbox.comparison.diff_models`. Use for cross-toolchain semantic-equality checks. |
| **Conditions** | [`yeastgem.conditions`](yeastgem/conditions.py) | `apply(model, name)` — minimal_Y6, anaerobic, glycine_nitrogen, nitrogen_limitation. Files under [`data/conditions/`](../../data/conditions/). |
| **Biomass** | [`yeastgem.biomass`](yeastgem/biomass.py) | `sum_biomass`, `scale_biomass`, `rescale_pseudoreaction`, `set_gam`, `change_amino_acid_ratio`. Configured from [`data/yeastgem/ids.yml`](../../data/yeastgem/ids.yml). |
| **Annotations** | [`yeastgem.missing_fields`](yeastgem/missing_fields.py) | `add_sbo_terms`, `load_delta_g`, `save_delta_g`. |
| **Model tests** | [`yeastgem.model_tests`](yeastgem/model_tests/) | `growth` (Tobias 2013 chemostat R²), `essential_genes` (Stanford KO collection), `anaerobic_flux_predictions`, `plot_anaerobic`, `find_duplicated_rxns`. |
| **Curation** | [`yeastgem.curation`](yeastgem/curation.py) | `curate_mets_rxns_genes` / `..._from_tsv` — batch curation from data tables with the yeast `s_`/`r_` id prefixes. |

## Layout

```
code/python/
  yeastgem/                # the package
    io.py                  # read_yeast_model + commit_yeast_model
    compare.py             # backwards-compat shim → raven_toolbox.comparison
    config.py              # YeastIDs loader (data/yeastgem/ids.yml)
    conditions.py          # apply(model, name)
    biomass.py             # sum_biomass / scale_biomass / set_gam / AA-ratio
    missing_fields.py      # add_sbo_terms, ΔG CSV persistence
    curation.py            # batch curation wrapper
    model_tests/           # Tier-3 benchmarks (growth, essential genes, …)
  tests/                   # pytest suite (65 tests across the package)
    reference/             # MATLAB-produced verification artefacts +
                           #   the runPhase*.m drivers
  pyproject.toml
  PORTING_PLAN.md
  UPSTREAM_CANDIDATES.md
```

## Running the tests

```bash
cd code/python
pytest -q
```

Tests load the real model once per session (~2 min) and exercise every
public function on it. ruff is the linter:

```bash
ruff check code/python
```

The CI workflow under
[`.github/workflows/python.yml`](../../.github/workflows/python.yml)
runs the same checks across Python 3.10 / 3.11 / 3.12, plus two
cross-language parity gates (level-1 SBML round-trip vs the committed
model, level-2 metric parity vs the committed reference values).

## Where work happens

Code under [`yeastgem/`](yeastgem/) is *only* yeast-specific
configuration and orchestration. Anything generic — model diff,
condition application, biomass scaling, curation, annotation — lives
in [raven-toolbox](https://github.com/SysBioChalmers/raven-toolbox).
Functions tracked for future upstreaming are in
[UPSTREAM_CANDIDATES.md](UPSTREAM_CANDIDATES.md).
