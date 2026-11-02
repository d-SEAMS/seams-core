import os
from pathlib import Path

import pytest

ROOT = Path(os.environ.get("DSEAMS_ROOT", Path(__file__).resolve().parents[2]))


@pytest.fixture(scope="session")
def root() -> Path:
    return ROOT


@pytest.fixture(scope="session")
def tip4p(root: Path) -> Path:
    return root / "input" / "traj" / "exampleTraj.lammpstrj"


@pytest.fixture(scope="session")
def cubic(root: Path) -> Path:
    return root / "input" / "traj" / "mW_cubic.lammpstrj"


@pytest.fixture(scope="session")
def si(root: Path) -> Path:
    return root / "input" / "traj" / "genice_sI.lammpstrj"
