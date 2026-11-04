"""Locate the yeast-GEM repository root (standard library only)."""
from __future__ import annotations

import os
from pathlib import Path

try:
    from dotenv import find_dotenv  # optional, kept for backwards compat
except ImportError:  # pragma: no cover - dotenv is a soft dependency
    find_dotenv = None  # type: ignore[assignment]


def find_repo_root() -> Path:
    """Locate the yeast-GEM repo root.

    Resolution order:

    1. The ``YEAST_GEM_PATH`` environment variable, if set.
    2. Walk up from this file looking for ``model/yeast-GEM.yml``.
    3. ``find_dotenv`` (historical convention; .env at repo root).
    4. Walk up from CWD looking for ``model/yeast-GEM.yml``.

    The marker is the ``.yml``, not ``.xml``: the yml is curators' source
    of truth and always tracked, while the xml is a generated artifact
    that may not exist at all on a fresh checkout (yeast-GEM#379 stage 2).
    """
    override = os.environ.get("YEAST_GEM_PATH")
    if override:
        return Path(override).resolve()
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "model" / "yeast-GEM.yml").exists():
            return parent
    if find_dotenv is not None:
        env = find_dotenv(usecwd=True)
        if env:
            return Path(env).parent.resolve()
    cwd = Path.cwd().resolve()
    for parent in (cwd, *cwd.parents):
        if (parent / "model" / "yeast-GEM.yml").exists():
            return parent
    raise FileNotFoundError(
        "Cannot locate the yeast-GEM repository root. "
        "Set the YEAST_GEM_PATH environment variable to the repo root, "
        "or place a .env file there."
    )
