<!-- https://docs.dseams.info/reference/julia-api.md -->

# Julia {#julia}

Documenter-style docstrings from [Dseams.jl](https://github.com/d-SEAMS/Dseams.jl) and the `Seams` package in that repository. `Dseams` is the language front. `Seams` calls `libseams` with `ccall`.

## Native API {#native-api}

Generated from the sources named below. A name on this page is a public call in that language.

### Dseams.Base {#dseams-base}

```julia
Base
```

One frame. `ptr` is owned by the finalizer, which calls `seams_frame_free`.

### Dseams.cages {#dseams-cages}

```julia
cages(frame; cutoff=3.5, type=0, k=4, graph="seeded") Hexagonal and double-diamond membership
```

cages(frame; cutoff=3.5, type=0, k=4, graph=”seeded”) Hexagonal and double-diamond membership. `cage` is 0 water, 1 hexagonal only, 2 double-diamond only, 3 both. The four counts are disjoint and match `seams cages`. `cutoff` is angstrom. A frame with no atoms gives an empty `cage` and four zero counts.

### Dseams.chill_plus {#dseams-chill-plus}

```julia
chill_plus(frame; cutoff = 3.5, type = 0, graph = "cutoff") CHILL+ labels
```

chill_plus(frame; cutoff = 3.5, type = 0, graph = “cutoff”) CHILL+ labels. `graph = “cutoff”` (the default) bonds every pair closer than `cutoff`. `graph = “mutual”` bonds the mutual four nearest neighbours; `cutoff + 1.5` only sizes the search cells and does not bound a bond, so an atom with no neighbour within `cutoff` can still be bonded. A frame with no atoms gives empty `labels` and `names` and nine zero counts.

### Dseams.engine_version {#dseams-engine-version}

```julia
engine_version() Version of the engine behind the loaded `libseams_lang`
```

engine_version() Version of the engine behind the loaded `libseams_lang`. ```jldoctest julia> using Dseams julia> VersionNumber(Dseams.engine_version()) >= v”2.14.0” true ```

### Dseams.example_path {#dseams-example-path}

```julia
example_path(name) Path of a trajectory shipped with the package
```

example_path(name) Path of a trajectory shipped with the package. `ice1h.lammpstrj` is a 128-molecule ice Ih cell, `exampleTraj.lammpstrj` a TIP4P/ice frame and `nucleation_f1.lammpstrj` a frame of an mW nucleation run.

### Dseams.features {#dseams-features}

```julia
features() Optional reader backends supported by the loaded `libseams_lang`, returned as `(chemfiles = Bool, readcon = Bool)`
```

features() Optional reader backends supported by the loaded `libseams_lang`, returned as `(chemfiles = Bool, readcon = Bool)`. ```jldoctest julia> using Dseams julia> keys(features()) == (:chemfiles, :readcon) true ```

### Dseams.has_chemfiles {#dseams-has-chemfiles}

```julia
has_chemfiles() `true` when the loaded engine was built with chemfiles support (PDB, GRO, DCD, CIF, and related trajectory formats), `false` when only the native LAMMPS and XYZ readers are built in
```

has_chemfiles() `true` when the loaded engine was built with chemfiles support (PDB, GRO, DCD, CIF, and related trajectory formats), `false` when only the native LAMMPS and XYZ readers are built in. ```jldoctest julia> using Dseams julia> has_chemfiles() isa Bool true ```

### Dseams.has_readcon {#dseams-has-readcon}

```julia
has_readcon() `true` when the loaded engine was built with `readcon-core` support for `.con` trajectories, `false` otherwise
```

has_readcon() `true` when the loaded engine was built with `readcon-core` support for `.con` trajectories, `false` otherwise. ```jldoctest julia> using Dseams julia> has_readcon() isa Bool true ```

### Dseams.ice_ih_path {#dseams-ice-ih-path}

```julia
ice_ih_path() Path of the Ice Ih example, `example_path("ice1h.lammpstrj")`.
```

ice_ih_path() Path of the Ice Ih example, `example_path(“ice1h.lammpstrj”)`.

### Dseams.jump_tau90 {#dseams-jump-tau90}

```julia
jump_tau90(first, second, dt=1.0; oxygen_type=1, hydrogen_type=2) Time for an H-H vector to rotate 90 degrees between two frames
```

jump_tau90(first, second, dt=1.0; oxygen_type=1, hydrogen_type=2) Time for an H-H vector to rotate 90 degrees between two frames. Returns `dt` when any water rotates that far, and `NaN` otherwise: when none does, when a frame has no water, or when `dt` is not positive.

### Dseams.library_path {#dseams-library-path}

```julia
library_path() The `libseams_lang` that Dseams opens: `ENV["SEAMS_LANG_LIB"]`, else the `libseams_lang` preference, else the one from `libseams_jll`
```

library_path() The `libseams_lang` that Dseams opens: `ENV\[“SEAMS_LANG_LIB”\]`, else the `libseams_lang` preference, else the one from `libseams_jll`. Where the JLL has no build, `”libseams_lang”` is looked up on the loader path.

### Dseams.module {#dseams-module}

```julia
module
```

Dseams Julia front for the d-SEAMS engine. Calls `libseams_lang` with `ccall`. The library comes from `libseams_jll`. A local build is used instead when the `SEAMS_LANG_LIB` environment variable or the `libseams_lang` preference names it; see \[`Dseams.set_library!`\](@ref). CxxWrap is not required. Ice Ih (`data/ice1h.lammpstrj`, every atom) is 128 hexagonal waters and 256 primitive six-rings. The C API names re-exported here call `libseams` (`SEAMS_LIB`, else the loader path) and are System, sea

### Dseams.mutable {#dseams-mutable}

```julia
mutable
```

Frame One frame from `libseams_lang`. The finalizer calls `seams_frame_free`. Do not use `ptr` after the frame is finalized.

### Dseams.neighbour_pairs {#dseams-neighbour-pairs}

```julia
neighbour_pairs(frame; cutoff=3.5, type=0, graph="cutoff") Directed pairs in the neighbour list
```

neighbour_pairs(frame; cutoff=3.5, type=0, graph=”cutoff”) Directed pairs in the neighbour list. The self entry is omitted. A frame with no atoms, or no pair within `cutoff`, has 0 pairs. `graph = “mutual”` is the four nearest neighbours and is not bounded by `cutoff`.

### Dseams.nframes {#dseams-nframes}

```julia
nframes(path) Number of frames in a trajectory
```

nframes(path) Number of frames in a trajectory. Native LAMMPS dumps and XYZ files are always supported; PDB, GRO, and other chemfiles formats are counted when \[`has_chemfiles`\](@ref) is `true`.

### Dseams.nop {#dseams-nop}

```julia
nop(frame) Number of atoms in `frame`.
```

nop(frame) Number of atoms in `frame`.

### Dseams.primitive_rings {#dseams-primitive-rings}

```julia
primitive_rings(frame; cutoff=3.5, type=0, graph="mutual", max_depth=6) Histogram of primitive rings
```

primitive_rings(frame; cutoff=3.5, type=0, graph=”mutual”, max_depth=6) Histogram of primitive rings. Index with the ring size: `rings\[6\]`, which is 0 when there is no such ring. `max_depth` is at least 3, else \[`SeamsError`\](@ref). A frame with no atoms or no ring gives an empty histogram.

### Dseams.q4_q6 {#dseams-q4-q6}

```julia
q4_q6(frame; cutoff=3.5, type=0) Local and Lechner-Dellago q4 and q6 on the cutoff shell
```

q4_q6(frame; cutoff=3.5, type=0) Local and Lechner-Dellago q4 and q6 on the cutoff shell. `cutoff` is angstrom. An atom with an empty shell has 0 for all four values, and a frame with no atoms gives empty vectors.

### Dseams.read_frame {#dseams-read-frame}

```julia
read_frame(path; frame=1, type=0) Read one frame
```

read_frame(path; frame=1, type=0) Read one frame. Native LAMMPS trajectory dumps and XYZ files are always supported; PDB, GRO, DCD, and CIF require \[`has_chemfiles`\](@ref), and `.con` files require \[`has_readcon`\](@ref). `frame` is 1-based. `type = 0` guesses oxygen. A positive type keeps that species. A negative type keeps every atom. A frame that holds no atoms, or none of the kept `type`, is a valid frame with `nop(frame) == 0`, and every analysis accepts it. A `frame` below 1, or past the la

### Dseams.set_library! {#dseams-set-library}

```julia
set_library!(path) set_library!(nothing) Use the `libseams_lang` at `path` (for instance `build/bindings` of a seams-core checkout) instead of the one in `libseams_jll`, or go back to the JLL with `nothing`
```

set_library!(path) set_library!(nothing) Use the `libseams_lang` at `path` (for instance `build/bindings` of a seams-core checkout) instead of the one in `libseams_jll`, or go back to the JLL with `nothing`. Stored in `LocalPreferences.toml` of the active project; restart Julia to take effect.

### Dseams.signature_cages {#dseams-signature-cages}

```julia
signature_cages(frame, spec; cutoff=3.5, type=0, graph="cutoff") Closed cages whose ring census matches `spec`, for example `"512"`
```

signature_cages(frame, spec; cutoff=3.5, type=0, graph=”cutoff”) Closed cages whose ring census matches `spec`, for example `”512”`. A frame with no atoms or no closed cage has 0 cages and 0 atoms.

### Dseams.steinhardt {#dseams-steinhardt}

```julia
steinhardt(frame, order_l; cutoff=3.5, type=0, graph=nothing) Local Steinhardt parameter
```

steinhardt(frame, order_l; cutoff=3.5, type=0, graph=nothing) Local Steinhardt parameter. `order_l` is 3, 4, 6, 8, or 12, and any other degree throws \[`SeamsError`\](@ref). Omitting `graph` uses the cutoff shell. `”mutual”` is the four nearest neighbours, not bounded by `cutoff`. An atom with an empty shell has `ql` and `qlbar` equal to 0, and a frame with no atoms gives empty vectors.

### Dseams.struct {#dseams-struct}

```julia
struct
```

CHILL+ labels, per-atom names, and the nine-label histogram.

### Dseams.symmetry {#dseams-symmetry}

```julia
symmetry(frame; shell_rmsd=0.12) Per-atom shell label
```

symmetry(frame; shell_rmsd=0.12) Per-atom shell label. 0 is ice Ih, 1 is ice Ic, 2 is stacking disorder, and -1 is anything else, which includes an atom with no neighbours. A frame with no atoms gives an empty vector.

### Seams.chill_plus {#seams-chill-plus}

```julia
chill_plus(sys; cutoff=3.5, graph=CUTOFF, rule=PUBLISHED) chill_plus(xyz, lengths; kwargs...) CHILL+ on a [`System`](@ref), or on a `3 × n` matrix in an orthorhombic cell
```

chill_plus(sys; cutoff=3.5, graph=CUTOFF, rule=PUBLISHED) chill_plus(xyz, lengths; kwargs…) CHILL+ on a \[`System`\](@ref), or on a `3 × n` matrix in an orthorhombic cell. `cutoff` is angstrom.

### Seams.const {#seams-const}

```julia
const
```

Published CHILL+ rule. Matches `SEAMS_CHILL_PUBLISHED`.

### Seams.module {#seams-module}

```julia
Seams Julia front for the d-SEAMS C API (`include/seams.h`)
```

Seams Julia front for the d-SEAMS C API (`include/seams.h`). Set `SEAMS_LIB` to `libseams.so` when it is not on the loader path. Handles are mutable structs with finalizers. A failing call throws \[`SeamsError`\](@ref). The message is `seams_handle_error`. `System`, `Neighbours`, `Rings`, and `CageSet` from the included C layer own their pointers. Result structs iterate in field order, so `labels, counts = chill_plus(sys)` still unpacks. Ice Ih (seams-core `input/traj/ice1h.lammpstrj`, every atom)

### Seams.num_threads {#seams-num-threads}

```julia
num_threads() Width the next analysis call will request.
```

num_threads() Width the next analysis call will request.

### Seams.q4q6 {#seams-q4q6}

```julia
q4q6(sys; cutoff=3.5) Local q4 and q6 on the cutoff shell
```

q4q6(sys; cutoff=3.5) Local q4 and q6 on the cutoff shell. `cutoff` is angstrom.

### Seams.set_num_threads {#seams-set-num-threads}

```julia
set_num_threads(n) OpenMP width for later analysis
```

set_num_threads(n) OpenMP width for later analysis. `n = 0` restores the saved default. A negative width throws \[`SeamsError`\](@ref). The width is process-global.

### Seams.set_ortho! {#seams-set-ortho}

```julia
set_ortho!(sys, xyz, lengths) Store a `3 × n` coordinate matrix in an orthorhombic cell at the origin
```

set_ortho!(sys, xyz, lengths) Store a `3 × n` coordinate matrix in an orthorhombic cell at the origin. `lengths` is `(lx, ly, lz)` in angstrom. The library copies the arrays.

### Seams.six_rings {#seams-six-rings}

```julia
six_rings(sys; cutoff=3.5) Primitive six-ring count and the fraction of atoms on a six-ring
```

six_rings(sys; cutoff=3.5) Primitive six-ring count and the fraction of atoms on a six-ring. The ring handle is freed before this function returns.

### Seams.struct {#seams-struct}

```julia
struct
```

CHILL+ labels and the nine population counts. Iterates as `(labels, counts)`.

### Seams.version {#seams-version}

```julia
version() Version text from `seams_version_string`
```

version() Version text from `seams_version_string`. The caller does not free it.

## C API {#c-api}

Each native call and the `seams.h` or `seams_lang.h` function it uses. The C column is the stable ABI.

| Native | C |
| --- | --- |
| Dseams.Base | no direct C call |
| Dseams.cages | seams_frame_cages |
| Dseams.chill_plus | seams_frame_chill_plus |
| Dseams.engine_version | seams_lang_version, seams_last_error |
| Dseams.example_path | no direct C call |
| Dseams.features | no direct C call |
| Dseams.has_chemfiles | seams_frame_free |
| Dseams.has_readcon | seams_frame_free |
| Dseams.ice_ih_path | no direct C call |
| Dseams.jump_tau90 | seams_frames_jump_tau90 |
| Dseams.library_path | no direct C call |
| Dseams.module | no direct C call |
| Dseams.mutable | seams_frame_free |
| Dseams.neighbour_pairs | seams_frame_neighbour_pairs, seams_frame_neighbours |
| Dseams.nframes | seams_path_nframes |
| Dseams.nop | seams_frame_chill_plus, seams_frame_nop |
| Dseams.primitive_rings | seams_frame_rings |
| Dseams.q4_q6 | seams_frame_q4_q6 |
| Dseams.read_frame | seams_frame_read |
| Dseams.set_library! | no direct C call |
| Dseams.signature_cages | seams_frame_signature |
| Dseams.steinhardt | seams_frame_steinhardt, seams_frame_steinhardt_graph |
| Dseams.struct | no direct C call |
| Dseams.symmetry | seams_frame_symmetry |
| Seams.chill_plus | seams_chill_plus, seams_system_count |
| Seams.const | seams_set_num_threads, seams_version_string |
| Seams.module | seams_set_num_threads, seams_version_string |
| Seams.num_threads | seams_get_num_threads, seams_system_set_ortho |
| Seams.q4q6 | seams_q4q6, seams_system_count |
| Seams.set_num_threads | seams_set_num_threads |
| Seams.set_ortho! | seams_system_set_ortho |
| Seams.six_rings | seams_rings_destroy, seams_rings_find, seams_rings_fraction, seams_rings_histogram |
| Seams.struct | no direct C call |
| Seams.version | seams_set_num_threads, seams_version_string |
