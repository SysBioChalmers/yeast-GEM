"""Yeast-GEM condition presets — yeast-specific wrapper over raven-python.

The generic mechanism (prelude, cofactor pseudoreaction edits,
biomass stoichiometry deltas, bounds, uptake-count check) lives in
:func:`raven_python.conditions.apply_condition`. yeast-GEM contributes:

1. The condition data files under ``data/conditions/``.
2. The :func:`load_condition` helper that resolves a name to a path.
3. The ``amino_acid_ratio`` pre-step that calls
   :func:`yeastgem.biomass.change_amino_acid_ratio` (used by the
   anaerobic condition).
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import cobra
from raven_python.conditions import (
    apply_condition as _ra_apply_condition,
)
from raven_python.conditions import (
    load_condition as _ra_load_condition,
)

from yeastgem.io import REPO_PATH

_CONDITIONS_DIR = REPO_PATH / "data" / "conditions"


def load_condition(name: str, *, conditions_dir: Path | None = None) -> dict[str, Any]:
    """Load ``data/conditions/<name>.yml`` as a plain dict."""
    base = conditions_dir or _CONDITIONS_DIR
    path = base / f"{name}.yml"
    return _ra_load_condition(path)


def apply(model: cobra.Model, name: str) -> cobra.Model:
    """Apply the named yeast-GEM condition to ``model`` in place.

    Pipeline:

    1. Resolve ``name`` to the YAML file under ``data/conditions/``.
    2. If the file declares ``amino_acid_ratio``, run the yeast-specific
       :func:`yeastgem.biomass.change_amino_acid_ratio` first (mirrors
       the MATLAB ``applyYeastCondition.m``).
    3. Hand the parsed config to
       :func:`raven_python.conditions.apply_condition` for the generic
       prelude / cofactor / biomass-delta / bounds steps.
    """
    cfg = load_condition(name)
    if "amino_acid_ratio" in cfg:
        # Imported lazily to keep yeastgem.conditions free of
        # biomass-module overhead when callers only use the bound-diff
        # path.
        from yeastgem.biomass import change_amino_acid_ratio

        aerobic = cfg["amino_acid_ratio"] == "aerobic"
        change_amino_acid_ratio(model, aerobic=aerobic)
    _ra_apply_condition(model, cfg)
    return model
