# Julia

`bindings/julia/Dseams` is the package CI runs. It calls `libseams_lang`. Point `SEAMS_LANG_LIB` at that shared library and `DSEAMS_ROOT` at the checkout. `bindings/julia/Seams` is the thin `libseams` example.

```julia
using Dseams
frame = read_frame(joinpath(ENV["DSEAMS_ROOT"], "input/traj/mW_cubic.lammpstrj"); type = 1)
ice = chill_plus(frame)
println(ice.counts["cubic"])
found = cages(frame)
println(found.n_ic)
```

`chill_plus`, `cages`, `primitive_rings`, `signature_cages`, `symmetry`, `q4_q6`, and `steinhardt` are the frame calls. `steinhardt` with no `graph` uses the cutoff shell.
