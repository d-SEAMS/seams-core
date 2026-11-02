# compute dseams

In-run CHILL+ labels from the `seams_chill_plus_graph` C ABI.

```
compute ice all dseams
dump 1 all custom 1 dump.chill id type x y z c_ice
```

The optional keywords are `graph cutoff|mutual` and `cutoff R`
(angstrom). The defaults are `graph cutoff cutoff 3.5`, the list
`seams chill-plus` uses. `graph mutual` keeps the mutual four nearest
neighbours inside `R + 1.5`, as `seams chill-plus --graph mutual` does.

The dump column is the same integer table as `seams chill-plus`
(0 cubic, 1 hexagonal, 2 water, 3 interfacial, 4 clathrate,
5 interClathrate, 6 unclassified). Copy these two sources into
`LAMMPS/src` and add `-I` to the seams-core include path.

`in.dseams` reads `water.data` (four oxygen sites) with `pair_style zero`
and writes `dump.chill`. The compute copies `nlocal` xyz and calls
`seams_chill_plus_graph`; it does not take the LAMMPS neighbour list.

`tests/test_phase.cpp` checks that the compute's label function is
`seams_chill_plus`. `build_compute.sh PREFIX` clones LAMMPS
`stable_22Jul2025_update5`, copies these two sources into `src/`,
links `libyodaLib` from `PREFIX`, and runs `in.dseams`.
