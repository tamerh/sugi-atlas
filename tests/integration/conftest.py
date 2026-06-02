"""Fixtures for the corpus-level integration suite (see _harness for the loader
and the run instructions / dist resolution)."""
import pytest

from ._harness import _maybe_skip, _load, DIST


@pytest.fixture(scope="session")
def pages():
    _maybe_skip()
    ps = _load()
    assert ps, "no pages loaded from the dist"
    return ps


@pytest.fixture(scope="session")
def dist_root():
    _maybe_skip()
    return DIST
