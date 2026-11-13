# Julia

`bindings/julia/Dseams` is the package CI runs. It calls `libseams_lang`. Point `SEAMS_LANG_LIB` at that shared library and `DSEAMS_ROOT` at the checkout. `bindings/julia/Seams` is the thin `libseams` example.

`meson test docs_api_julia` runs `docs/examples/julia_api.jl` with
`SEAMS_LANG_LIB` pointed at `libseams_lang`.
`docs/examples/minimal_julia.jl` reads the
GenIce sI cell and calls `steinhardt`, `chill_plus`, and `primitive_rings`.
`meson test docs_examples` runs it.

```{literalinclude} ../../../docs/examples/minimal_julia.jl
:language: julia
```

```{literalinclude} ../../../docs/examples/julia_api.jl
:language: julia
```

`Frame`, `nop`, `read_frame`, `nframes`, `chill_plus`, `cages`, `primitive_rings`, `signature_cages`, `symmetry`, `q4_q6`, `jump_tau90`, and `steinhardt` are the exported `Dseams` calls. `steinhardt` with no `graph` uses the cutoff shell. `Seams` exports `System`, `version`, `set_num_threads`, `num_threads`, `set_ortho!`, `chill_plus`, `q4q6`, and `six_rings` over `libseams`.

`bindings/julia/Dseams/test/runtests.jl` and `bindings/julia/Seams/test/runtests.jl` are the tests `language_fronts` runs.

```{literalinclude} ../../../bindings/julia/Dseams/test/runtests.jl
:language: julia
:lines: 14-37,57-67
```

```{literalinclude} ../../../bindings/julia/Seams/test/runtests.jl
:language: julia
:lines: 19-36
```
