# Julia tutorial

`bindings/julia/Dseams` is the package CI tests. It calls `libseams_lang` and covers rings, cages, and the TIP4P counts. `bindings/julia/Seams` is the thin `libseams` example (`examples/ice.jl` builds an eight-atom diamond cell).

The cubic mW dump is ice Ic. CHILL+ labels all 4096 oxygens cubic.

```julia
using Dseams
root = get(ENV, "DSEAMS_ROOT", pwd())
frame = read_frame(joinpath(root, "input", "traj", "mW_cubic.lammpstrj"); type = 1)
ice = chill_plus(frame)
println(ice.counts["cubic"])
```

```bash
export DSEAMS_ROOT=$PWD
export SEAMS_LANG_LIB=$PWD/build/bindings/libseams_lang.so
export LD_LIBRARY_PATH=$PWD/build/bindings:$PWD/build/src
julia --project=bindings/julia/Dseams -e 'using Pkg; Pkg.instantiate()'
julia --project=bindings/julia/Dseams bindings/julia/Dseams/test/runtests.jl
```

The suite checks the mixed TIP4P frame, seeded cages, and that cubic count of 4096.
