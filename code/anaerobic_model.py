"""Alias of :mod:`apply_anaerobic`, so existing ``from anaerobic_model import
anaerobic_model`` imports keep working. New code should use
``from apply_anaerobic import apply_anaerobic``.
"""
from apply_anaerobic import anaerobic_model, apply_anaerobic  # noqa: F401
