import pytest
import os

@pytest.fixture(scope='session')
def files_path() -> str:
    return os.path.join(os.path.dirname(__file__), 'files')
