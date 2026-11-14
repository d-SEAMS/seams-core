# Dseams

Julia front for the [d-SEAMS](https://docs.dseams.info) engine: CHILL+,
primitive rings, ice cages, Steinhardt parameters and shell symmetry on
LAMMPS trajectories. The engine comes from `libseams_jll`.

Dseams and `libseams_jll` are not in the General registry yet. From a
seams-core checkout, develop the stand-in JLL and name a local build of
`libseams_lang`:

```sh
julia --project=bindings/julia/Dseams bindings/julia/Dseams/dev/develop_jll.jl
export SEAMS_LANG_LIB=$PWD/build/bindings/libseams_lang.so
```

```julia
using Dseams
f = read_frame(Dseams.ice_ih_path(); type = -1)
chill_plus(f; graph = "cutoff").counts["hexagonal"]   # 128
```

`docs/` has the manual (`julia --project=docs docs/make.jl`), including
how to use a local engine build instead of the JLL and how Dseams relates
to the `Seams` example next to it.

## Registering

`Project.toml` is ready for Registrator (`uuid`, `version`, and
`[compat]`, with `libseams_jll` held at `=2.14.0`). General can resolve
that dependency only after `libseams_jll` is registered from the
Yggdrasil recipe in `packaging/yggdrasil`. Until then, `dev/develop_jll.jl`
develops the stand-in JLL or a local `--deploy=local` build.
