"""yeastgem — Python interface to the yeast-GEM consensus model.

Python-side counterpart to the MATLAB functions under
[code/](../). See [PORTING_PLAN.md](../PORTING_PLAN.md) for scope and
[UPSTREAM_CANDIDATES.md](../UPSTREAM_CANDIDATES.md) for what may move
upstream later.
"""
from __future__ import annotations

from yeastgem import biomass, conditions, curation, model_tests
from yeastgem.compare import ComparisonReport, compare_models
from yeastgem.config import YeastIDs, load_ids
from yeastgem.io import (
    MODEL_PATH,
    REPO_PATH,
    commit_yeast_model,
    read_yeast_model,
    write_yeast_model,
)
from yeastgem.missing_fields import add_sbo_terms, load_delta_g, save_delta_g

__all__ = [
    "MODEL_PATH",
    "REPO_PATH",
    "ComparisonReport",
    "YeastIDs",
    "add_sbo_terms",
    "biomass",
    "commit_yeast_model",
    "compare_models",
    "conditions",
    "curation",
    "load_delta_g",
    "load_ids",
    "model_tests",
    "read_yeast_model",
    "save_delta_g",
    "write_yeast_model",
]

__version__ = "0.0.1.dev0"
