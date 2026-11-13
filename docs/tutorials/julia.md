# Julia tutorial

`bindings/julia/Dseams` is the package CI tests. It calls `libseams_lang` and covers rings, cages, and the TIP4P counts. `bindings/julia/Seams` is the thin `libseams` example (`examples/ice.jl` builds an eight-atom diamond cell).

The cubic mW dump is ice Ic. CHILL+ labels all 4096 oxygens cubic.

```{literalinclude} ../../../bindings/julia/Dseams/test/runtests.jl
:language: julia
:lines: 57-62
```

```bash
export DSEAMS_ROOT=$PWD
export SEAMS_LANG_LIB=$PWD/build/bindings/libseams_lang.so
export LD_LIBRARY_PATH=$PWD/build/bindings:$PWD/build/src
julia --project=bindings/julia/Dseams -e 'using Pkg; Pkg.instantiate()'
julia --project=bindings/julia/Dseams bindings/julia/Dseams/test/runtests.jl
```

The suite checks the mixed TIP4P frame, seeded cages, and that cubic count of 4096.


## Structure still

```{figure} ../tutorials/molecules/julia/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug julia` writes the 2400 by 1800 PNG.
```

