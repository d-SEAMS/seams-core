# LAMMPS

Meson builds `run_hook`, the in-run hook. It links `liblammps` and registers a
`fix external` callback. Every step the callback calls the same C API as
`seams chill-plus`. `engines/lammps/compute_dseams.cpp` is a per-atom style
you can copy into a LAMMPS source tree. This build does not compile it.

```
fix ext all external pf/callback 2 2
run 8
```

`in.dseams` is the deck. `ice.data` is 64 oxygens of cubic ice. Forces are
zeroed so the cell does not melt. The hook prints `overhead_s`, the mean
seconds spent inside d-SEAMS per step. `run_hook.cpp` includes `library.h`,
so the LAMMPS include directory is on the compile line.
