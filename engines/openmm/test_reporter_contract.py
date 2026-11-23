#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""The reporter descriptor keeps the tuple contract accepted by OpenMM 8.0."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "adapters" / "python"))

from dseams_insitu.openmm import DseamsReporter  # noqa: E402


class Simulation:
    currentStep = 5


reporter = DseamsReporter(object(), interval=4)
descriptor = reporter.describeNextReport(Simulation())
if descriptor != (3, True, False, False, False):
    raise SystemExit("OpenMM 8.0 reporter descriptor is {!r}".format(descriptor))
print("OpenMM 8.0 reporter tuple", descriptor)
