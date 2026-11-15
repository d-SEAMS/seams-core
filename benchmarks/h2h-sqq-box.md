# SQQ and the triclinic sH cell

SQQ 0.5.8 reads only an orthorhombic cell. The raw log is benchmarks/h2h-sqq-box.txt. With one reader configuration, `inspect_lammps_dump` accepts `input/traj/genice_sI.lammpstrj` and rejects `input/traj/genice_sH.lammpstrj`:

    Triclinic/tilted LAMMPS dump boxes are not supported: genice_sH.lammpstrj frame 0.

That dump's box header is `ITEM: BOX BOUNDS xy xz yz pp pp pp` with a nonzero xy tilt. The same call with `pbc.box_mode` set to `triclinic` is rejected before any frame is read:

    pbc.box_mode must be orthorhombic.

The checks live in the library, not in this harness.

`sqq/io/input/lammps.py` opens with "Strict orthorhombic LAMMPS trajectory input for SQQ." `inspect_lammps_dump` raises `ValueError` when the box header contains `xy`, `xz`, `yz`, `abc`, `avec`, `bvec` or `cvec`. `validate_lammps_data_box` raises the same class of error for a general triclinic DATA box and for a nonzero xy/xz/yz tilt. `sqq/io/input/trajectory.py` `parse_gro_box` rejects a 9-value GRO box whose tilt entries are nonzero and tells the caller to convert the frame to an orthorhombic representation first.

`sqq/config/scopes.py` requires `pbc.box_mode` to be `orthorhombic`. `sqq/config/serialize.py` documents that key as "currently orthorhombic only". `sqq/core/common/pbc.py` `minimum_image` is the orthorhombic minimum image, using the first three box entries. `sqq/core/sqq_py/spatial.py` is the cutoff search for orthorhombic coordinates.

The sH comparison in benchmarks/h2h-cages.md is therefore d-SEAMS alone. SQQ's cage counts are reported for the orthorhombic sI and sII cells.
