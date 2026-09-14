"""Shared fixtures for the bhmiepy test suite."""

from pathlib import Path

import pytest

from bhmiepy._upstream_cli import find_ref_binary

REPO_ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = REPO_ROOT / "upstream" / "examples"


@pytest.fixture(scope="session")
def bhmie_ref():
    """The pristine upstream CLI built from the upstream/ submodule.

    Fails loudly (no silent skip) when the submodule is uninitialized or
    the binary has not been built.
    """
    if not EXAMPLES.is_dir():
        raise RuntimeError(
            "upstream/ submodule is not initialized; run: "
            "git submodule update --init --recursive"
        )
    return find_ref_binary(REPO_ROOT)
