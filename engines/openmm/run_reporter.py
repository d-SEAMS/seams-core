# Copyright (c) d-SEAMS developers
"""Run the OpenMM reporter, or skip when OpenMM is not installed."""

import importlib.util
import runpy
import sys
from pathlib import Path

if importlib.util.find_spec("openmm") is None:
    print("openmm is not installed")
    sys.exit(77)

here = Path(__file__).resolve().parent
sys.path.insert(0, str(here))
runpy.run_path(str(here / "test_reporter.py"), run_name="__main__")
