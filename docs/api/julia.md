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

`chill_plus`, `cages`, `primitive_rings`, `signature_cages`, `symmetry`, `q4_q6`, `jump_tau90`, `nframes`, and `steinhardt` are the `Dseams` calls. `steinhardt` with no `graph` uses the cutoff shell. `Seams` exports `System`, `version`, `set_num_threads`, `num_threads`, `set_ortho!`, `chill_plus`, `q4q6`, and `six_rings` over `libseams`.

`bindings/julia/Dseams/test/runtests.jl` and `bindings/julia/Seams/test/runtests.jl` are the tests `language_fronts` runs.

```{literalinclude} ../../../bindings/julia/Dseams/test/runtests.jl
:language: julia
:lines: 14-37,57-67
```

```{literalinclude} ../../../bindings/julia/Seams/test/runtests.jl
:language: julia
:lines: 19-36
```
