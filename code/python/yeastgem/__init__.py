"""yeastgem — Python interface to the yeast-GEM consensus model.

Python-side counterpart to the MATLAB functions under
[code/](../). See [PORTING_PLAN.md](../PORTING_PLAN.md) for scope and
[UPSTREAM_CANDIDATES.md](../UPSTREAM_CANDIDATES.md) for what may move
upstream later.
"""
from __future__ import annotations

import importlib

# Submodules and names are imported when first used, so that
# ``import yeastgem.conditions`` does not import raven-toolbox (``yeastgem.io``
# and the modules built on it do).
_SUBMODULES = (
    "annotate", "biomass", "compare", "conditions", "config", "curation", "io",
    "missing_fields", "model_tests", "paths",
)
_ATTRIBUTES = {
    "annotate_gem": "annotate",
    "derive_annotation_tsvs": "annotate",
    "ComparisonReport": "compare",
    "compare_models": "compare",
    "YeastIDs": "config",
    "load_ids": "config",
    "MODEL_PATH": "io",
    "REPO_PATH": "io",
    "YAML_PATH": "io",
    "commit_yeast_model": "io",
    "load_yeast_yaml": "io",
    "read_yeast_model": "io",
    "save_yeast_yaml": "io",
    "write_yeast_model": "io",
    "add_sbo_terms": "missing_fields",
    "load_delta_g": "missing_fields",
    "save_delta_g": "missing_fields",
}


def __getattr__(name: str):
    if name in _SUBMODULES:
        return importlib.import_module(f"yeastgem.{name}")
    if name in _ATTRIBUTES:
        value = getattr(importlib.import_module(f"yeastgem.{_ATTRIBUTES[name]}"), name)
        globals()[name] = value
        return value
    raise AttributeError(f"module 'yeastgem' has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted({*globals(), *_SUBMODULES, *_ATTRIBUTES})


__all__ = [
    "MODEL_PATH",
    "REPO_PATH",
    "YAML_PATH",
    "ComparisonReport",
    "YeastIDs",
    "add_sbo_terms",
    "annotate",
    "annotate_gem",
    "biomass",
    "commit_yeast_model",
    "compare_models",
    "conditions",
    "curation",
    "derive_annotation_tsvs",
    "load_delta_g",
    "load_ids",
    "load_yeast_yaml",
    "model_tests",
    "read_yeast_model",
    "save_delta_g",
    "save_yeast_yaml",
    "write_yeast_model",
]

__version__ = "0.0.1.dev0"
