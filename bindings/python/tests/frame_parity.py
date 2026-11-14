# Copyright (c) d-SEAMS developers
"""Frame results against pydseamslib 2.11.1 on the documented dumps.

The oracle is a separate interpreter. Set PYDSEAMS211 to that tree and
PYDSEAMS211_PYTHON to its interpreter. The comparison never puts that
tree on this process's import path.
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

ROOT = Path(os.environ.get("DSEAMS_SRC", Path(__file__).resolve().parents[3]))
TRAJ = ROOT / "input" / "traj"


def _ours():
    os.environ.pop("DSEAMS_DOC_BUILD", None)
    os.environ["PYTHONPATH"] = str(ROOT / "bindings" / "python") + os.pathsep + os.environ.get(
        "PYTHONPATH", ""
    )
    sys.path.insert(0, str(ROOT / "bindings" / "python"))
    import pydseams as ds
    from pydseams.band import ring_census

    return ds, ring_census


def _oracle():
    root = os.environ.get("PYDSEAMS211", "")
    exe = os.environ.get("PYDSEAMS211_PYTHON", "")
    if not root or not exe:
        pytest.skip("PYDSEAMS211 and PYDSEAMS211_PYTHON name the 2.11.1 tree")
    script = r"""
import json, os, tempfile
from pathlib import Path
import ase.io
import pydseams as ds
from pydseams.band import ring_census
root = Path(os.environ["TRAJ"])
out = {}
frame = ds.read(str(root / "exampleTraj.lammpstrj"), atom_type=2)
out["chill"] = dict(frame.chill_plus())
cages = frame.cages()
out["cages"] = {"n_ih": cages.n_ih, "n_ic": cages.n_ic, "n_water": cages.n_water}
out["rings_default"] = ring_census(frame, 8)
stein = {}
for ell in (3, 4, 6):
    sample = frame.steinhardt(ell)
    stein[str(ell)] = {"ql": sample["ql"], "ql_bar": sample["ql_bar"]}
out["steinhardt"] = stein
td = tempfile.mkdtemp()
radius, values = frame.rdf_2d(td, cutoff=6.0, binwidth=0.5)
out["rdf"] = {"r": radius, "g": values}
out["mono_default"] = {str(k): v for k, v in frame.monolayer_rings(td, 1.0, max_depth=6).items()}
out["prisms"] = frame.find_prisms(output_dir=td, max_depth=4)
guest = ds.read(str(root / "genice_sI.lammpstrj"), atom_type=1)
found = guest.cages_by_signature("512")
verts = set()
for cage in found:
    verts.update(int(v) for v in cage.get("vertices", []))
out["sig"] = {"n": len(found), "verts": len(verts)}
mw = ds.read(str(root / "mW_cubic.lammpstrj"), frame=1, atom_type=1)
out["mw"] = {"cubic": int(mw.chill_plus().cubic), "n_ic": int(mw.cages().n_ic), "n": int(mw.n_atoms)}
cut = ds.read(str(root / "exampleTraj.lammpstrj"), atom_type=2, bonded="cutoff")
out["rings_cutoff"] = ring_census(cut, 8)
td2 = tempfile.mkdtemp()
out["mono_cutoff"] = {str(k): v for k, v in cut.monolayer_rings(td2, 1.0, max_depth=6).items()}
atoms = ase.io.read(str(root / "exampleTraj.lammpstrj"), format="lammps-dump-text")
atoms.set_chemical_symbols(["O" if int(t) == 2 else "H" for t in atoms.arrays["type"]])
ase_frame = ds.from_ase(atoms)
out["ase_n"] = int(ase_frame.n_atoms)
out["ase_chill"] = dict(ase_frame.chill_plus())
labelled = ase_frame.to_ase()
out["ase_symbols"] = int((labelled.symbols == "O").sum())
system = frame.to_solvis()
out["solvis"] = type(system).__name__
print(json.dumps(out))
"""
    env = os.environ.copy()
    env["PYTHONPATH"] = root
    env["TRAJ"] = str(TRAJ)
    env.pop("DSEAMS_DOC_BUILD", None)
    done = subprocess.run([exe, "-c", script], env=env, text=True, capture_output=True)
    if done.returncode != 0:
        raise AssertionError(done.stderr[-2000:] or done.stdout[-2000:])
    return json.loads(done.stdout.strip().splitlines()[-1])


def test_frame_matches_pydseamslib_211():
    ds, ring_census = _ours()
    oracle = _oracle()
    frame = ds.read(str(TRAJ / "exampleTraj.lammpstrj"), atom_type=2)
    ice = frame.chill_plus()
    # The documented cutoff census is 9 interClathrate and 241 water.
    # pydseamslib 2.11.1 chill_plus() has no graph argument and uses its
    # default neighbour list, which is the older 12 and 238.
    assert ice.interClathrate == 9
    assert ice.water == 241
    assert oracle["chill"]["interClathrate"] == 12
    assert oracle["chill"]["water"] == 238
    cages = frame.cages()
    assert (cages.n_ih, cages.n_ic, cages.n_water) == (
        oracle["cages"]["n_ih"],
        oracle["cages"]["n_ic"],
        oracle["cages"]["n_water"],
    )
    assert cages.n_water == 250
    # Default 2.11.1 ring_census follows bonded="hbond" whenever the file
    # still names hydrogens. The page describes the cutoff neighbour list.
    # Forcing bonded="cutoff" on 2.11.1 reproduces seams_frame_rings.
    census = ring_census(frame, max_size=8)
    assert census == {int(k): v for k, v in oracle["rings_cutoff"].items()}
    assert census == {3: 37, 4: 273, 5: 34, 6: 167}
    assert {int(k): v for k, v in oracle["rings_default"].items()} == {3: 1, 4: 230, 5: 4, 6: 192}
    for ell in (3, 4, 6):
        sample = frame.steinhardt(ell)
        other = oracle["steinhardt"][str(ell)]
        assert sample["ql"] == pytest.approx(other["ql"])
        assert sample["ql_bar"] == pytest.approx(other["ql_bar"])
    with tempfile.TemporaryDirectory() as tmp:
        radius, values = frame.rdf_2d(tmp, cutoff=6.0, binwidth=0.5)
        assert radius == pytest.approx(oracle["rdf"]["r"])
        assert values == pytest.approx(oracle["rdf"]["g"], rel=1e-4, abs=1e-4)
        mono = frame.monolayer_rings(tmp, sheet_area=1.0, max_depth=6)
        cutoff = {int(k): v for k, v in oracle["mono_cutoff"].items()}
        for size, row in cutoff.items():
            assert mono[size]["count"] == row["count"]
            assert mono[size]["coverage_xy"] == pytest.approx(row["coverage_xy"], rel=1e-4, abs=1e-3)
        # The default 2.11.1 sheet uses the hydrogen-bond rings, so its
        # counts are the 1/230/4/192 census, not the cutoff kernel.
        default = {int(k): v["count"] for k, v in oracle["mono_default"].items()}
        assert default[3] == 1 and default[4] == 230 and default[6] == 192
        rows = frame.find_prisms(output_dir=tmp, max_depth=4)
    assert [row["size"] for row in rows] == [3, 4]
    assert oracle["prisms"] is None
    guest = ds.read(str(TRAJ / "genice_sI.lammpstrj"), atom_type=1)
    found = guest.cages_by_signature("512")
    assert len(found) == oracle["sig"]["n"] == 2
    # seams_frame_signature returns the cage count and the distinct
    # vertex count. It does not return one vertex list per cage.
    assert found.n_atoms == oracle["sig"]["verts"] == 40
    ice = ds.read(str(TRAJ / "mW_cubic.lammpstrj"), frame=1, atom_type=1)
    assert ice.n_atoms == oracle["mw"]["n"] == 4096
    assert ice.chill_plus().cubic == oracle["mw"]["cubic"] == 4096
    assert ice.cages().n_ic == oracle["mw"]["n_ic"] == 4096
    import ase.io

    atoms = ase.io.read(str(TRAJ / "exampleTraj.lammpstrj"), format="lammps-dump-text")
    atoms.set_chemical_symbols(["O" if int(t) == 2 else "H" for t in atoms.arrays["type"]])
    ase_frame = ds.from_ase(atoms)
    assert ase_frame.n_atoms == oracle["ase_n"] == 250
    ase_ice = ase_frame.chill_plus()
    # from_ase keeps the selected oxygens and classifies that cutoff frame.
    # That is the published 9 and 241. pydseamslib 2.11.1 still builds its
    # hydrogen-bond list from the hydrogens left in the Atoms object, so
    # its from_ase census stays on the older 12 and 238 neighbour list.
    assert ase_ice.interClathrate == 9
    assert ase_ice.water == 241
    assert oracle["ase_chill"]["interClathrate"] == 12
    assert oracle["ase_chill"]["water"] == 238
    labelled = ase_frame.to_ase()
    assert int((labelled.symbols == "O").sum()) == oracle["ase_symbols"]
    assert "ice_type" in labelled.arrays
    system = frame.to_solvis()
    assert type(system).__name__ == oracle["solvis"]
