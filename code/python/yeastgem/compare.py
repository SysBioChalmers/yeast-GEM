"""Backwards-compatibility shim for the cross-language model comparator.

The implementation moved to ``raven_python.comparison.diff`` (the
``diff_models`` function) as part of the phase-3.5 restructure. This
module re-exports it under the original ``compare_models`` /
``ComparisonReport`` names so existing yeast-GEM callers keep working;
new code should import from ``raven_python.comparison`` directly.
"""
from __future__ import annotations

from raven_python.comparison import (
    DEFAULT_ANNOTATION_KEYS,
)
from raven_python.comparison import (
    DiffReport as ComparisonReport,
)
from raven_python.comparison import (
    diff_models as compare_models,
)

__all__ = [
    "DEFAULT_ANNOTATION_KEYS",
    "ComparisonReport",
    "compare_models",
]
