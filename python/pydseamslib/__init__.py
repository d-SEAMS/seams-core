"""d-SEAMS Python front end.

Load a frame and ask for ice::

    import pydseamslib as ds
    frame = ds.read("water.lammpstrj")
    print(frame.chill_plus())
    print(frame.cages())

ASE Atoms work the same way::

    frame = ds.from_ase(atoms)          # default: oxygen
    atoms = frame.to_ase()
"""

from . import _core
from .frame import CageScore, Frame, IceCounts, read

__version__ = "2.0.0"

# Drop-in name used in the 2.0 docs and tests
Trajectory = Frame


def from_ase(atoms, select="O", cutoff=3.5, bonded="auto"):
    """Build a Frame from an ASE Atoms. select is a symbol or atomic number."""
    return Frame.from_ase(atoms, select=select, cutoff=cutoff, bonded=bonded)


def from_arrays(positions, cell, numbers=None, cutoff=3.5, bonded="cutoff"):
    """Build a Frame from (N, 3) positions and three box lengths."""
    return Frame.from_arrays(
        positions, cell, numbers=numbers, cutoff=cutoff, bonded=bonded
    )


__all__ = [
    "Frame",
    "Trajectory",
    "IceCounts",
    "CageScore",
    "read",
    "from_ase",
    "from_arrays",
    "_core",
    "__version__",
]
