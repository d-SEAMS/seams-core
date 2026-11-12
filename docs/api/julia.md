# Julia

`bindings/julia/Seams` calls `libseams` with `ccall`. Point `SEAMS_LIB` at the
shared library.

```julia
using Seams
sys = System()
set_ortho!(sys, xyz, (lx, ly, lz))  # xyz is 3 by n
labels, counts = chill_plus(sys; rule=Seams.PUBLISHED)
q4, q6 = q4q6(sys)
six, fraction = six_rings(sys)
```

`set_num_threads`, `num_threads`, and `version` match the C API.
`examples/ice.jl` prints the cubic count.
