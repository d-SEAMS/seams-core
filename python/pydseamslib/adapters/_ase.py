"""Deprecated. Use pydseamslib.from_ase / Frame.to_ase."""

from pydseamslib import from_ase as _from_ase


def to_frame(atoms, select="O", cutoff=3.5, bonded="auto"):
    return _from_ase(atoms, select=select, cutoff=cutoff, bonded=bonded)
