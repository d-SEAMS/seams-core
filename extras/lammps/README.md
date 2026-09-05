# compute dseams

In-run CHILL+ labels from the `seams_chill_plus` C ABI.

```
compute ice all dseams
dump 1 all custom 1 dump.chill id type x y z c_ice
```

The dump column is the same integer table as `seams chill-plus`
(0 cubic, 1 hexagonal, 2 water, 3 interfacial, 4 clathrate,
5 interClathrate, 6 unclassified). Copy these two sources into
`LAMMPS/src` and add `-I` to the seams-core include path.

`tests/test_phase.cpp` checks that the compute's label function is
`seams_chill_plus`.
