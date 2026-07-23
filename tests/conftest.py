"""Shared fixtures for the pipeline test suite.

The pipeline scripts live under scripts/ and are named with a leading digit
(e.g. 1_filter.py), so they can't be imported with a normal `import` statement.
`import_script` loads them by file path instead.
"""
import importlib.util
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = PROJECT_ROOT / "scripts"


def import_script(filename):
    """Import a pipeline script (e.g. '1_filter.py') as a module by file path."""
    module_name = filename[:-3]
    spec = importlib.util.spec_from_file_location(module_name, SCRIPTS_DIR / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def filter_module():
    return import_script("1_filter.py")


@pytest.fixture(scope="session")
def nearby_module():
    return import_script("3_find_nearby_mutations.py")


@pytest.fixture(scope="session")
def pipeline_utils_module():
    """pipeline_utils.py isn't digit-prefixed, so unlike the fixtures above it
    doesn't need import_script's path-loader workaround -- a plain import
    after adding scripts/ to sys.path is exactly what every production
    script (1_filter.py, 3_find_nearby_mutations.py, etc.) already does.
    """
    sys.path.insert(0, str(SCRIPTS_DIR))
    import pipeline_utils
    return pipeline_utils


@pytest.fixture(scope="session")
def download_module():
    return import_script("2_download_structures.py")


class FakeResponse:
    """Minimal stand-in for requests.Response used to mock UniProt API calls."""

    def __init__(self, text, headers=None, status_code=200):
        self.text = text
        self.headers = headers or {}
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")
