"""Shared pytest fixtures for yeastgem tests."""
from __future__ import annotations

import pytest

from yeastgem import read_yeast_model


@pytest.fixture(scope="session")
def model():
    """The yeast-GEM model loaded once per test session."""
    return read_yeast_model()
