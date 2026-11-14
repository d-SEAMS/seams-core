# Dseams

Julia front for the [d-SEAMS](https://docs.dseams.info) engine: CHILL+,
primitive rings, ice cages, Steinhardt parameters and shell symmetry on
LAMMPS trajectories. The engine comes from `libseams_jll`.

```julia
using Dseams
f = read_frame(Dseams.ice_ih_path(); type = -1)
chill_plus(f; graph = "cutoff").counts["hexagonal"]   # 128
```

`docs/` has the manual (`julia --project=docs docs/make.jl`), including
how to use a local engine build instead of the JLL and how Dseams relates
to the `Seams` example next to it.
