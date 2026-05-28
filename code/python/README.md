# yeastgem — Python port of the yeast-GEM functions

This directory hosts the in-development Python interface to yeast-GEM. It
is the Python counterpart to the MATLAB code under [../](..).

## Status

Early scaffolding. See [PORTING_PLAN.md](PORTING_PLAN.md) for the full
plan and [UPSTREAM_CANDIDATES.md](UPSTREAM_CANDIDATES.md) for the
function-level upstream tracking.

## Install (development)

From a yeast-GEM checkout:

```bash
pip install -e code/python/[dev]
```

This installs the `yeastgem` package (cobrapy-based; no ravengem
dependency by design).

## Quick start

```python
from yeastgem import read_yeast_model

model = read_yeast_model()           # cobra.Model
print(model.optimize().objective_value)
```

`yeastgem` auto-detects the repo root via the package location, the
`YEAST_GEM_PATH` environment variable, or a `.env` file at the repo
root (historical convention) — in that order.

## Layout

```
code/python/
  yeastgem/          # the package
    io.py            # read/write the model (commit_yeast_model lands later)
  tests/             # pytest unit tests
  pyproject.toml
  PORTING_PLAN.md
  UPSTREAM_CANDIDATES.md
```
