<!-- https://docs.dseams.info/fronts/lammps.md -->

# LAMMPS tutorial {#lammps-tutorial}

`engines/lammps/in.dseams` is a frozen cubic cell of 64 atoms. `fix ext` uses `pf/callback 2 2`, so the library calls d-SEAMS every two steps. Forces are cleared and the cell stays cubic.

```default
fix ext all external pf/callback 2 2
thermo 2
timestep 0.001
run 8
```

The driver is `run_hook`, built as `build/engines/run_hook`. It registers the callback before the `run` line and checks that successive calls are one stride apart. Positions stay in angstroms.

```bash
meson compile -C build run_hook
cd engines/lammps
LD_LIBRARY_PATH=../../build/src ../../build/engines/run_hook in.dseams
```

A passing run prints a line that starts `lammps calls 5 stride 2 ice 320`. `overhead_s` on that line is the mean seconds per call. Five visits of 64 atoms are all cubic. The next line is `lammps hexagonal 0 six-rings 128 q6 0.628539`.

`compute dseams` in `engines/lammps/compute_dseams.cpp` is the same CHILL+ call for a LAMMPS build that compiles styles in.
