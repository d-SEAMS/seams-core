# Copyright (c) d-SEAMS developers
"""Checks on an installed pydseamslib: version, compatibility imports, libseams_lang.

Not named test_* so the in-tree unittest run, which has no installed
libseams_lang, does not collect it. The wheel and conda tests pass it
to pytest by path.
"""

import math

import numpy as np

import dseams


def diamond(nrep=3):
    a = 2.75 * 4.0 / math.sqrt(3.0)
    fcc = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5))
    basis = ((0, 0, 0), (0.25, 0.25, 0.25))
    points = [
        ((i + f[0] + b[0]) * a, (j + f[1] + b[1]) * a, (k + f[2] + b[2]) * a)
        for i in range(nrep)
        for j in range(nrep)
        for k in range(nrep)
        for f in fcc
        for b in basis
    ]
    length = nrep * a
    return np.asarray(points, dtype=np.float64), (length, length, length)


def test_version():
    assert dseams.__version__ == dseams.version()
    parts = dseams.__version__.split(".")
    assert len(parts) == 3 and all(part.isdigit() for part in parts)


def test_compat_imports():
    import pydseams
    import pydseamslib

    assert pydseamslib.__version__ == dseams.__version__
    assert pydseamslib.Frame is dseams.Frame
    assert callable(pydseams.read)


def test_lang_library_resolves():
    from dseams import userframe

    path = userframe._lang_library()
    assert path.endswith((".so", ".dylib", ".dll")), path


def test_cubic_ice_shell_and_cages():
    pos, box = diamond()
    labels = dseams.symmetry(pos, box)
    assert labels.shape == (len(pos),)
    # Every site of a perfect diamond lattice is ice Ic.
    assert np.all(labels == 1), np.unique(labels, return_counts=True)
    per, (n_hex, n_cubic, n_water, n_both) = dseams.cages(pos, box)
    assert per.shape == (len(pos),)
    assert n_cubic == len(pos)
    assert n_hex == 0
