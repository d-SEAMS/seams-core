# LAMMPS tutorial

`engines/lammps/in.dseams` is a real LAMMPS run. `fix ext` uses `pf/callback 2 2`, so the library calls d-SEAMS every two steps. Forces are cleared and the cell stays cubic.

```{literalinclude} ../../../engines/lammps/in.dseams
:lines: 12-15
```

The driver is `run_hook`. It registers the callback before the `run` line and checks that successive calls are one stride apart.

```bash
meson compile -C build run_hook
cd engines/lammps
LD_LIBRARY_PATH=../../build/src ../../build/engines/run_hook in.dseams
```

A passing run prints `lammps calls 5 stride 2 ice 320`. Five visits of 64 atoms are all cubic. `compute dseams` in `compute_dseams.cpp` is the same call for a LAMMPS build that compiles styles in.
