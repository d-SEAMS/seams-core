"""Bond-orientational order parameters from pydseamslib.

``Frame.steinhardt`` calls ``yoda.steinhardtQl`` on the cutoff
neighbour list. The CLI ``seams steinhardt`` command uses the mutual
four-nearest list, so the means on the mixed frame are not the CLI
means.
"""

import math

import pytest

pydseams = pytest.importorskip("pydseams")


def _mean(values):
    finite = [float(value) for value in values if math.isfinite(float(value))]
    assert finite
    return sum(finite) / len(finite)


def test_q3_means_on_the_cutoff_list(tip4p):
    frame = pydseams.read(tip4p, atom_type=2)
    q3 = frame.steinhardt(3)
    assert _mean(q3["ql"]) == pytest.approx(0.374897, abs=1e-4)
    assert _mean(q3["ql_bar"]) == pytest.approx(0.187542, abs=1e-4)
    q12 = frame.steinhardt(12)
    assert _mean(q12["ql"]) == pytest.approx(0.512940, abs=1e-4)


def test_q6_is_finite_on_cubic_ice(cubic):
    frame = pydseams.read(cubic, atom_type=1)
    q6 = frame.steinhardt(6)
    assert len(q6["ql"]) == 4096
    assert math.isfinite(_mean(q6["ql"]))
