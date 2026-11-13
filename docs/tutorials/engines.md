# Engines tutorial

Eight molecules of cubic ice are `engines/tutorial/ice.data`. Lengths are angstroms. `meson test engines_tutorial` reads that file and classifies it with `libseams`.

LAMMPS `run_hook` takes the coordinates LAMMPS already stores in angstroms and calls `seams_chill_plus`. The deck next to the data file is the same hook the in-tree test uses, pointed at this cell:

```{literalinclude} ../../../engines/tutorial/in.dseams
```

The OpenMM plugin reads nanometres. It multiplies positions by 10 before the same call. The tutorial test builds those nanometres from the data file, scales them back, and requires the same eight cubic labels.

```{literalinclude} ../../../engines/tutorial_walk.py
:language: python
:pyobject: openmm_angstroms
```

`engines/openmm/test_reporter.py` is the full reporter when OpenMM is installed. `meson test lammps_hook` and `meson test openmm_reporter` run when those libraries were found at configure time. Both end at `seams_chill_plus` on a cubic cell.

```bash
meson test -C build engines_tutorial --print-errorlogs
```

A passing run prints `engines tutorial ok atoms 8 cubic 8`.
