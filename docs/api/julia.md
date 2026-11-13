# Julia

`bindings/julia/Dseams` calls `libseams_lang`. `bindings/julia/Seams` calls
`libseams` (`include/seams.h`). Set `SEAMS_LANG_LIB` and `SEAMS_LIB`.
`DSEAMS_ROOT` is the checkout.

Both modules are packages with `Project.toml`, exports, and `Pkg.test`.
Handles are mutable structs with finalizers. Result values are structs
with `Base.show`. A failing call throws `SeamsError`. The message is
`seams_handle_error` for the C API and `seams_last_error` for `libseams_lang`.

## Quickstart

```julia
using Dseams
frame = read_frame(path; type = 2)
ice = chill_plus(frame)
println(ice)
println(ice.counts["water"])
```

On `input/traj/exampleTraj.lammpstrj` with oxygen type 2 the cutoff shell
is 9 interClathrate and 241 water. The mutual shell is 1 clathrate,
7 interClathrate, and 242 water. Mean q3 on the cutoff shell is 0.374897.
`primitive_rings` indexes by size: `rings[6] == 224` at depth 8.

`Seams` on the cubic cell:

```julia
using Seams
set_num_threads(1)
sys = System()
set_ortho!(sys, xyz, (edge, edge, edge))
labels, counts = chill_plus(sys; rule = PUBLISHED)
q4, q6 = q4q6(sys)
six, fraction = six_rings(sys)
```

`ChillLabels`, `Q4Q6Pair`, and `SixRings` iterate in field order, so the
unpacking above is the struct itself. Every atom is cubic. Mean q6 is
above 0.3. The six-ring fraction is above 0.5.

## Ownership

`Frame` calls `seams_frame_free` from its finalizer. `System` calls
`seams_system_destroy`. A null pointer is ignored, so the finalizer is
safe to run twice. Neighbour and ring handles created inside `six_rings`,
`seams_cages`, and `seams_rings_vertices` are destroyed before those
functions return. `show(sys)` prints the atom count.

## Errors

`SeamsError` is the exception. Pass a path that the reader rejects and
the message is the handle string. An empty `System` has `seams_chill`
throw before any label buffer is trusted.

## Multiple dispatch

`chill_plus` accepts a `Frame`, a trajectory path, or, in `Seams`, a
`System` and a `3 × n` matrix. `read_frame`, `nop`, and `nframes` stay
on the path and frame methods.

## Native reference

Documenter builds the reference from the docstrings:

```bash
julia --project=bindings/julia/Dseams/docs -e 'using Pkg; Pkg.instantiate()'
julia --project=bindings/julia/Dseams/docs bindings/julia/Dseams/docs/make.jl
julia --project=bindings/julia/Seams/docs -e 'using Pkg; Pkg.instantiate()'
julia --project=bindings/julia/Seams/docs bindings/julia/Seams/docs/make.jl
```

That is the same split readcon uses between the bindings page and the
generated crate pages. `doctest = false` in `make.jl` because the
numerical checks live in `Pkg.test`, which needs the shared libraries.

```{literalinclude} ../../../docs/examples/minimal_julia.jl
:language: julia
```

```{literalinclude} ../../../docs/examples/julia_api.jl
:language: julia
```

```{literalinclude} ../../../bindings/julia/Dseams/test/runtests.jl
:language: julia
:lines: 14-37,57-67
```

```{literalinclude} ../../../bindings/julia/Seams/test/runtests.jl
:language: julia
:lines: 19-36
```

<!-- API_SYMBOLS_BEGIN -->

## Symbol reference

### Frame

Purpose: One frame read from a dump or built in memory.

Parameters: none

Units: object

Example: `docs/examples/minimal_julia.jl`

### cages

Purpose: Hexagonal and double-diamond cage counts.

Parameters: frame; cutoff

Units: cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_julia.jl`

### chill_plus

Purpose: CHILL+ label of each atom.

Parameters: frame; cutoff

Units: cutoff in angstrom; labels are dimensionless codes

Example: `docs/examples/minimal_julia.jl`

### jump_tau90

Purpose: Time for an H-H vector to rotate 90 degrees between two frames.

Parameters: first; second; dt

Units: dt and the result in the caller's time unit

Example: `docs/examples/minimal_julia.jl`

### nframes

Purpose: Number of frames in a trajectory path.

Parameters: path

Units: dimensionless count

Example: `docs/examples/minimal_julia.jl`

### nop

Purpose: Number of atoms in a frame.

Parameters: frame

Units: dimensionless count

Example: `docs/examples/minimal_julia.jl`

### primitive_rings

Purpose: Histogram of primitive rings.

Parameters: frame; cutoff; max_ring

Units: cutoff in angstrom; bin k is a dimensionless ring count

Example: `docs/examples/minimal_julia.jl`

### q4_q6

Purpose: Local q4 and q6.

Parameters: frame; cutoff

Units: cutoff in angstrom; q4 and q6 are dimensionless

Example: `docs/examples/minimal_julia.jl`

### read_frame

Purpose: Read one frame from a dump.

Parameters: path; frame; type

Units: frame index is dimensionless; coordinates are angstrom

Example: `docs/examples/minimal_julia.jl`

### signature_cages

Purpose: Closed cages matching a ring-census spec.

Parameters: frame; cutoff; spec

Units: cutoff in angstrom; spec is text

Example: `docs/examples/minimal_julia.jl`

### steinhardt

Purpose: Mean local Steinhardt q_l.

Parameters: frame; cutoff; order_l

Units: cutoff in angstrom; ql is dimensionless

Example: `docs/examples/minimal_julia.jl`

### symmetry

Purpose: Shell label for Ih, Ic, or stacking disorder.

Parameters: frame

Units: labels are dimensionless codes

Example: `docs/examples/minimal_julia.jl`
### System

Purpose: Call into libseams through ccall.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/julia_api.jl`

### seams_system_count

Purpose: Call into libseams through ccall.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/julia_api.jl`

### seams_system_read

Purpose: Call into libseams through ccall.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/julia_api.jl`

### seams_chill

Purpose: Call into libseams through ccall.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/julia_api.jl`

### seams_steinhardt

Purpose: Call into libseams through ccall.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/julia_api.jl`

### seams_cages

Purpose: Call into libseams through ccall.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/julia_api.jl`

### seams_signature_cages

Purpose: Call into libseams through ccall.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/julia_api.jl`

### seams_neighbours_knn

Purpose: Call into libseams through ccall.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/julia_api.jl`

### seams_rdf

Purpose: Call into libseams through ccall.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/julia_api.jl`

### seams_running_cn

Purpose: Call into libseams through ccall.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/julia_api.jl`

### seams_coordination

Purpose: Call into libseams through ccall.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/julia_api.jl`

### seams_hbonds

Purpose: Call into libseams through ccall.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/julia_api.jl`

### seams_density_axis

Purpose: Call into libseams through ccall.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/julia_api.jl`

### seams_fingerprint

Purpose: Call into libseams through ccall.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/julia_api.jl`

### seams_contact_pairs

Purpose: Call into libseams through ccall.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/julia_api.jl`

### seams_domains

Purpose: Call into libseams through ccall.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/julia_api.jl`

### seams_prisms

Purpose: Call into libseams through ccall.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/julia_api.jl`

### seams_monolayer

Purpose: Call into libseams through ccall.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/julia_api.jl`

### seams_rings_vertices

Purpose: Call into libseams through ccall.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/julia_api.jl`

### seams_site_table

Purpose: Call into libseams through ccall.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/julia_api.jl`

<!-- API_SYMBOLS_END -->
