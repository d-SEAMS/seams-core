# LAMMPS

`engines/lammps/compute_dseams.cpp` is a per-atom compute. It reads local
coordinates and the LAMMPS cell, calls `seams_chill_plus`, and stores the
label in a per-atom vector. Build it into LAMMPS, or compile it as an object
against the LAMMPS headers:

```
compute ice all dseams cutoff 3.5 graph cutoff rule published
```

The Ubuntu package used here does not include the PLUGIN package, so the
compute cannot be dlopened. The end-to-end test is `run_hook`, which links
`liblammps` and registers a `fix external` callback. The callback is the
in-run hook: every step it calls the same C API the compute calls.

`in.dseams` is the deck. `ice.data` is 64 oxygens of cubic ice. Forces are
zeroed so the cell does not melt. The hook prints `overhead_s`, the mean
seconds spent inside d-SEAMS per step.
