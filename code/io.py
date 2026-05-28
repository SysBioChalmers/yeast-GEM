"""DEPRECATED — this module moved to the ``yeastgem`` package.

``code/io.py`` is kept as a forwarding shim for backwards
compatibility. Install the new package with::

    pip install -e code/python/

and import from ``yeastgem`` instead::

    from yeastgem import read_yeast_model, write_yeast_model

This shim will be removed in a future release.
"""
import warnings

warnings.warn(
    "code/io.py is deprecated. Install with `pip install -e code/python/` "
    "and import from `yeastgem` instead.",
    DeprecationWarning,
    stacklevel=2,
)

from yeastgem.io import (  # noqa: E402, F401
    MODEL_PATH,
    REPO_PATH,
    read_yeast_model,
    write_yeast_model,
)
