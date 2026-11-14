# Dseams

Julia front for the [d-SEAMS](https://docs.dseams.info) engine: CHILL+,
primitive rings, ice cages, Steinhardt parameters and shell symmetry on
LAMMPS trajectories.

```julia
using Dseams
f = read_frame(Dseams.ice_ih_path(); type = -1)
chill_plus(f; graph = "cutoff").counts["hexagonal"]       # 128
primitive_rings(f; graph = "cutoff", max_depth = 6)[6]    # 256
```

## Dseams and Seams

seams-core has two Julia packages. **Dseams** is the one to install. It
reads trajectories and has the same frame analyses as the R package. It
calls `libseams_lang`, the C ABI shared by the R and Julia fronts, and is
what the bindings CI and the documentation run.

`bindings/julia/Seams` is a short example of the stable C API in
`include/seams.h` (`libseams`). It works on coordinates the caller passes
in and reads no files. It is kept as an example of that API and is not
registered.

## The engine library

`libseams_jll` provides the engine for Linux (glibc and musl; x86_64,
aarch64, powerpc64le, riscv64) and macOS. To use a local build of
seams-core instead, name its `libseams_lang`:

```julia
Dseams.set_library!("/path/to/seams-core/build/bindings/libseams_lang.so")
# restart Julia; Dseams.set_library!(nothing) goes back to the JLL
```

This writes the `libseams_lang` preference to `LocalPreferences.toml`.
The `SEAMS_LANG_LIB` environment variable takes precedence over the
preference. With either one set, `libseams_jll` is not loaded.

To replace the whole JLL artifact with an install prefix instead, use
an `Overrides.toml` in `~/.julia/artifacts`, as for any JLL.

## Example data

`Dseams.example_path(name)` returns one of the trajectories in `data/`:
`ice1h.lammpstrj` (ice Ih, 128 oxygens), `exampleTraj.lammpstrj` (a
TIP4P/ice frame) and `nucleation_f1.lammpstrj` (one frame of an mW
nucleation run). The `.source.txt` files there record where they came
from.
