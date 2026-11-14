# Julia tutorial

`bindings/julia/Dseams` is the package CI tests. It calls `libseams_lang` and covers rings, cages, and the TIP4P counts. `bindings/julia/Seams` is the thin `libseams` example (`examples/ice.jl` builds an eight-atom diamond cell).

The cubic mW dump is ice Ic. CHILL+ labels all 4096 oxygens cubic. `docs/examples/core_counts.jl` checks Ice Ih and sII. `meson test docs_core_counts` runs it.

```{literalinclude} ../../../docs/examples/core_counts.jl
:language: julia
```

```{literalinclude} ../../../bindings/julia/Dseams/test/runtests.jl
:language: julia
:lines: 82-86
```

```bash
export DSEAMS_ROOT=$PWD
export SEAMS_LANG_LIB=$PWD/build/bindings/libseams_lang.so
export LD_LIBRARY_PATH=$PWD/build/bindings:$PWD/build/src
julia --project=bindings/julia/Dseams bindings/julia/Dseams/dev/develop_jll.jl
julia --project=bindings/julia/Dseams bindings/julia/Dseams/test/runtests.jl
```

The suite checks the mixed TIP4P frame, seeded cages, and that cubic count of 4096.
The cubic dump is in `input/traj`, so that block runs only with `DSEAMS_ROOT` set.
`Dseams` and `libseams_jll` are not in the General registry yet, so `Pkg.instantiate()` alone cannot resolve the JLL.
`dev/develop_jll.jl` develops the stand-in in `bindings/julia/Dseams/dev/libseams_jll`, which has no binaries, and instantiates; `SEAMS_LANG_LIB` names the engine.


## Structure still

```{figure} ../tutorials/molecules/julia/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug julia` writes the 2400 by 1800 PNG.
```

