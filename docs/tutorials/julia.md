# Julia tutorial

`bindings/julia/Seams` calls `libseams` with `ccall`. Julia 1.11 needs `Test` listed under `[extras]`.

```julia
using Seams
a = 2.75 * 4 / sqrt(3)
xyz = zeros(3, 8)
sys = System()
set_ortho!(sys, xyz, (a, a, a))
labels, counts = chill_plus(sys)
println(version(), " ", counts[1])
```

```bash
export SEAMS_LIB=$PWD/build/src/libseams.so
export LD_LIBRARY_PATH=$PWD/build/src
julia --project=bindings/julia/Seams -e 'using Pkg; Pkg.test()'
julia --project=bindings/julia/Seams bindings/julia/Seams/examples/ice.jl
```

The suite is eight checks on a 64-atom cubic cell: both CHILL+ rules, q6, and the six-ring fraction.
