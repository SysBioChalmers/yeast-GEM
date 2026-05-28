"""yeastgem — Python interface to the yeast-GEM consensus model.

Python-side counterpart to the MATLAB functions under
[code/](../). See [PORTING_PLAN.md](../PORTING_PLAN.md) for scope and
[UPSTREAM_CANDIDATES.md](../UPSTREAM_CANDIDATES.md) for what may move
upstream later.
"""
from __future__ import annotations

from yeastgem.compare import ComparisonReport, compare_models
from yeastgem.io import (
    MODEL_PATH,
    REPO_PATH,
    read_yeast_model,
    write_yeast_model,
)

__all__ = [
    "MODEL_PATH",
    "REPO_PATH",
    "ComparisonReport",
    "compare_models",
    "read_yeast_model",
    "write_yeast_model",
]

__version__ = "0.0.1.dev0"
