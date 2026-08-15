"""ASE Atoms in and out of a Frame. Optional: pip install 'pydseamslib[ase]'."""

from __future__ import annotations


def _require_ase():
    try:
        import ase
        from ase import Atoms
        from ase.data import chemical_symbols
    except ImportError as exc:
        raise ImportError(
            "ASE interop needs the ase package. "
            "Install it with: pip install 'pydseamslib[ase]'"
        ) from exc
    return ase, Atoms, chemical_symbols


def _mask(atoms, select):
    if select is None:
        return [True] * len(atoms)
    if isinstance(select, str):
        return [sym == select for sym in atoms.get_chemical_symbols()]
    number = int(select)
    return [int(z) == number for z in atoms.numbers]


def _cell_lengths(atoms):
    cell = atoms.get_cell()
    return [float(cell[0, 0]), float(cell[1, 1]), float(cell[2, 2])]


def frame_from_ase(cls, atoms, select="O", cutoff=3.5, bonded="auto"):
    _, Atoms, _ = _require_ase()
    if not hasattr(atoms, "get_positions"):
        raise TypeError("from_ase expects an ASE Atoms object")
    keep = _mask(atoms, select)
    if not any(keep):
        raise ValueError(
            f"no atoms matched select={select!r}; "
            f"symbols={sorted(set(atoms.get_chemical_symbols()))}"
        )
    positions = [
        xyz for xyz, yes in zip(atoms.get_positions(), keep) if yes
    ]
    numbers = [
        int(z) for z, yes in zip(atoms.numbers, keep) if yes
    ]
    symbols = [
        s for s, yes in zip(atoms.get_chemical_symbols(), keep) if yes
    ]
    origin = atoms.get_celldisp().reshape(3)
    from .frame import _cloud_from_positions

    cloud = _cloud_from_positions(
        positions, _cell_lengths(atoms), numbers, box_low=origin
    )
    h_pos = [
        xyz
        for xyz, s in zip(atoms.get_positions(), atoms.get_chemical_symbols())
        if s == "H"
    ]
    h_cloud = None
    if h_pos:
        h_cloud = _cloud_from_positions(
            h_pos,
            _cell_lengths(atoms),
            [1] * len(h_pos),
            box_low=origin,
        )
    if bonded == "auto":
        bonded = "hbond" if h_cloud is not None else "cutoff"
    return cls(
        atom_type=int(numbers[0]),
        cutoff=cutoff,
        bonded=bonded,
        cloud=cloud,
        h_cloud=h_cloud,
        symbols=symbols,
    )


def frame_to_ase(frame):
    _, Atoms, _ = _require_ase()
    n = frame.n_atoms
    if frame._symbols is not None:
        symbols = list(frame._symbols)
    else:
        # LAMMPS type 8 is oxygen in the periodic table; type 1 or 2
        # from a water dump is not. Prefer O for the analysed species.
        fallback = "O"
        symbols = [fallback] * n
    atoms = Atoms(
        symbols=symbols,
        positions=frame.positions,
        cell=frame.box,
        pbc=True,
    )
    ice = [pt.iceType.name for pt in frame.cloud.pts]
    if any(name != "unclassified" for name in ice):
        atoms.arrays["ice_type"] = ice
    atoms.arrays["hc"] = [False] * n
    atoms.info["dseams_n_atoms"] = n
    return atoms
