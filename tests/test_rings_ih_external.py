#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Ring counts on ideal ice Ih, where other codes use the same rings.

seams_rings_find returns each Franzblau primitive ring once. On this
crystal that is two hexagons per molecule and no smaller ring. It then
labels each hexagon as an HC basal, an HC prismatic, or a hexagon that
is neither. Those labels are not extra rings, and OVITO, cycless, and
rings-statistics-matter do not emit them.

The strict prism-block label is also absent: bulk ice Ih has no pair of
same-size rings with a perfect bond matching and an empty intersection.
The hexagonal cages are the HC basal and HC prismatic labels.

cycless ``cycles_iter`` with fractional coordinates, and
rings-statistics-matter ``ring_statistics`` with ``supercell`` on,
count the same hexagons here. The coordinate filter drops a cycle
that winds the cell. Without it, and with ``supercell`` left off on
this 128-molecule cell, both tools also count squares and extra
hexagons. Those are not primitive rings of the periodic graph.

The oxygen network has girth 6, so a hexagon has no chord, and this
lattice does not contain the spur that makes a chordless hexagon fail
the shortest-path test. That spur is in ``markdown/rings-theory.md``.
A chordless search and this search diverge there, not on ideal ice Ih.

OVITO ``FindRingsModifier`` runs on the 4096-molecule cell; on the
128-molecule cell the short axis lets it count image-winding squares
and extra hexagons. OVITO 3.17.2 reports all 8192 hexagons there.
OVITO 3.15.5, the locked linux-64 wheel, reports 6815 of them, most
away from the cell faces. So the check is that OVITO finds only
hexagons, each one once, each one of the hexagons cycless finds on
that cell, and from OVITO 3.17 on all 8192 of them. The comparison
runs when ``ovito`` imports. A missing package exits 77.
"""

from __future__ import annotations

import importlib.metadata
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

try:
    import numpy as np
    from gen_ideal_ice import ice_ih
except ImportError as exc:
    print(f"ice Ih generator skipped: {exc}", file=sys.stderr)
    raise SystemExit(77) from exc

CUTOFF = 3.5
MOLECULES = 128


def _positions():
    frame = ice_ih(MOLECULES)
    pos = np.asarray(frame["positions"], dtype=np.float64)
    lengths = np.diag(np.asarray(frame["lattice"], dtype=np.float64))
    if pos.shape != (MOLECULES, 3):
        raise SystemExit(f"ice Ih cell has {pos.shape[0]} molecules")
    return pos, lengths


def _graph(pos, lengths):
    import networkx as nx

    cut2 = CUTOFF * CUTOFF
    graph = nx.Graph()
    graph.add_nodes_from(range(len(pos)))
    for i in range(len(pos)):
        delta = pos[i + 1:] - pos[i]
        delta -= lengths * np.round(delta / lengths)
        dist2 = np.einsum("ij,ij->i", delta, delta)
        for j in np.nonzero((dist2 > 1e-8) & (dist2 <= cut2))[0]:
            graph.add_edge(i, i + 1 + int(j))
    return graph


def _expect(census, name):
    want = [0, 0, 0, 2 * MOLECULES]
    if list(census) != want:
        raise SystemExit(f"{name} census {census} is not {want}")
    try:
        version = importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        version = ""
    print(f"{name} {version} sizes 3-6 {want}")


def cycless_rings(pos, lengths):
    from cycless.cycles import cycles_iter

    # pos drops a cycle whose minimum-image steps wind the cell. Without
    # that filter the same graph also yields squares and extra hexagons.
    fractional = pos / lengths
    return [tuple(sorted(ring)) for ring in cycles_iter(_graph(pos, lengths), 6, pos=fractional)]


def cycless_census(pos, lengths):
    census = [0, 0, 0, 0]
    for ring in cycless_rings(pos, lengths):
        size = len(ring)
        if 3 <= size <= 6:
            census[size - 3] += 1
    return census


def matter_census(pos, lengths):
    from ase import Atoms
    from rings_statistics_matter import ring_statistics

    atoms = Atoms(
        symbols=["O"] * len(pos),
        positions=np.ascontiguousarray(pos, dtype=np.float64),
        cell=lengths.tolist(),
        pbc=True,
    )
    # The 128-molecule cell is narrower than the search depth. supercell
    # repeats it so a ring cannot meet its image, then maps the count back.
    counts, _rings = ring_statistics(
        atoms,
        cutoff=CUTOFF,
        maxlvl=8,
        supercell=True,
        deduplicate=True,
        backend="rust",
        num_threads=1,
    )
    values = np.asarray(counts, dtype=np.float64)
    return [int(round(float(v))) for v in values[2:6]]


def _ovito_rings(pos, lengths):
    """OVITO's census for sizes 3-6 and its hexagons as sorted index tuples."""
    import os

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    os.environ.setdefault("OVITO_THREAD_COUNT", "1")
    from ovito.data import DataCollection
    from ovito.modifiers import CreateBondsModifier, FindRingsModifier
    from ovito.pipeline import Pipeline, StaticSource

    mat = np.zeros((3, 4))
    mat[:, :3] = np.diag(lengths)
    data = DataCollection()
    if hasattr(data, "create_cell"):
        data.create_cell(mat, pbc=(True, True, True))
        parts = data.create_particles(count=len(pos))
    else:
        from ovito.data import Particles, SimulationCell

        cell = SimulationCell(pbc=(True, True, True))
        cell.matrix = mat
        data.objects.append(cell)
        parts = Particles()
        data.objects.append(parts)
    parts.create_property("Position", data=np.asarray(pos))
    pipe = Pipeline(source=StaticSource(data=data))
    pipe.modifiers.append(CreateBondsModifier(cutoff=CUTOFF))
    pipe.modifiers.append(FindRingsModifier(maximum_ring_size=6))
    out = pipe.compute()
    hist = out.tables["ring-size-histogram"]
    found = {int(s): int(c) for s, c in zip(hist["Ring Size"][...], hist["Count"][...])}
    hexagons = []
    if "6-rings" in out.tables.keys():
        members = np.asarray(out.tables["6-rings"]["Particle Indices"][...])
        hexagons = [tuple(sorted(int(i) for i in row)) for row in members]
    return [found.get(size, 0) for size in range(3, 7)], hexagons


def _ensure_gl():
    # OVITO's PySide6 looks up libGL by name. A system copy outside the
    # default loader path is enough; the process has to start again so
    # the loader sees it.
    # An interpreter that already loads OVITO keeps its own loader path.
    try:
        import ovito  # noqa: F401

        return
    except ImportError:
        pass
    gl = "/usr/lib/x86_64-linux-gnu"
    lib = os.path.join(gl, "libGL.so.1")
    if not os.path.isfile(lib) or os.environ.get("SEAMS_OVITO_LIBS") == "1":
        return
    current = os.environ.get("LD_LIBRARY_PATH", "")
    os.environ["LD_LIBRARY_PATH"] = gl + ((":" + current) if current else "")
    os.environ["SEAMS_OVITO_LIBS"] = "1"
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    os.execv(sys.executable, [sys.executable, *sys.argv])


def main():
    _ensure_gl()
    pos, lengths = _positions()
    ran = 0
    try:
        census = cycless_census(pos, lengths)
    except ImportError as exc:
        print(f"cycless skipped: {exc}", file=sys.stderr)
    else:
        _expect(census, "cycless")
        ran += 1
    try:
        census = matter_census(pos, lengths)
    except ImportError as exc:
        print(f"rings-statistics-matter skipped: {exc}", file=sys.stderr)
    else:
        _expect(census, "rings-statistics-matter")
        ran += 1
    try:
        # 4096 molecules. On the 128-molecule cell the short axis lets
        # FindRings count image-winding squares and extra hexagons.
        frame = ice_ih(4096)
        big = np.asarray(frame["positions"], dtype=np.float64)
        big_lengths = np.diag(np.asarray(frame["lattice"], dtype=np.float64))
        census, hexagons = _ovito_rings(big, big_lengths)
        import ovito
    except ImportError as exc:
        print(f"ovito skipped: {exc}", file=sys.stderr)
    else:
        if census[:3] != [0, 0, 0] or census[3] != len(hexagons):
            raise SystemExit(f"ovito census {census} has a ring that is not a listed hexagon")
        if len(set(hexagons)) != len(hexagons):
            raise SystemExit("ovito lists a hexagon twice")
        if tuple(ovito.version[:2]) >= (3, 17) and len(hexagons) != 2 * 4096:
            raise SystemExit(f"ovito {ovito.version} finds {len(hexagons)} hexagons, not {2 * 4096}")
        try:
            primitive = set(cycless_rings(big, big_lengths))
        except ImportError as exc:
            print(f"ovito ring membership skipped: {exc}", file=sys.stderr)
            primitive = None
        if primitive is not None:
            if len(primitive) != 2 * 4096:
                raise SystemExit(f"cycless finds {len(primitive)} hexagons at 4096, not {2 * 4096}")
            stray = [ring for ring in hexagons if ring not in primitive]
            if stray:
                raise SystemExit(f"ovito hexagon {stray[0]} is not a primitive ring")
        print(
            f"ovito {importlib.metadata.version('ovito')} hexagons {len(hexagons)} "
            f"of {2 * 4096}" + (", all primitive" if primitive is not None else "")
        )
        ran += 1
    if ran == 0:
        raise SystemExit(77)
    print(f"compared {ran} tools")


if __name__ == "__main__":
    main()
