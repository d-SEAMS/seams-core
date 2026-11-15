# Julia

`bindings/julia/Dseams` calls `libseams_lang`. `bindings/julia/Seams` calls
`libseams` (`include/seams.h`). Set `SEAMS_LANG_LIB` and `SEAMS_LIB`.
`DSEAMS_ROOT` is the checkout.

`Dseams` and its engine package `libseams_jll` are not in the General registry yet, so `Pkg.add("Dseams")` does not work. `bindings/julia/Dseams/dev/develop_jll.jl` develops a stand-in `libseams_jll` with no binaries, and `SEAMS_LANG_LIB` or `Dseams.set_library!` names a local `libseams_lang`. `packaging/yggdrasil/L/libseams` is the Yggdrasil recipe for the JLL.

`meson test docs_api_julia` runs `docs/examples/julia_api.jl` with
`SEAMS_LANG_LIB` pointed at `libseams_lang`.
`docs/examples/minimal_julia.jl` reads the
GenIce sI cell and calls `steinhardt`, `chill_plus`, and `primitive_rings`.
`meson test docs_examples` runs it.

Both modules are packages with `Project.toml`, exports, and `Pkg.test`.
Handles are mutable structs with finalizers. Result values are structs
with `Base.show`. A failing call throws `SeamsError`. The message is
`seams_handle_error` for the C API and `seams_last_error` for `libseams_lang`.

## Quickstart

```julia
using Dseams
path = joinpath(ENV["DSEAMS_ROOT"], "input", "traj", "exampleTraj.lammpstrj")
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
a = 2.75 * 4 / sqrt(3)
fcc = ((0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))
basis = ((0.0, 0.0, 0.0), (0.25, 0.25, 0.25))
sites = [((i, j, k) .+ f .+ b) .* a for i in 0:1, j in 0:1, k in 0:1, f in fcc, b in basis]
xyz = reduce(hcat, collect.(vec(sites)))
edge = 2a
set_num_threads(1)
sys = System()
set_ortho!(sys, xyz, (edge, edge, edge))
labels, counts = chill_plus(sys; rule = Seams.PUBLISHED)
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
The built reference is linked from this book:

```{raw} html
<p><a href="../julia-dseams/index.html">Dseams Documenter reference</a>
and <a href="../julia-seams/index.html">Seams Documenter reference</a>.</p>
```

`docs/examples/core_counts.jl` checks Ice Ih (128 hexagonal, 256
six-rings) and sII (136 clathrate, 144 five-rings, 16 six-rings).
`meson test docs_core_counts` runs it.

```{literalinclude} ../../../docs/examples/core_counts.jl
:language: julia
```

```{literalinclude} ../../../docs/examples/minimal_julia.jl
:language: julia
```

```{literalinclude} ../../../docs/examples/julia_api.jl
:language: julia
```

```{literalinclude} ../../../bindings/julia/Dseams/test/runtests.jl
:language: julia
:lines: 27-50,68-71,82-86
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

Returns: A Frame.

Ownership: The finalizer calls seams_frame_free.

Errors: Error when the pointer is null.

Units: object

Example: `docs/examples/minimal_julia.jl`

### Seams.CUTOFF

Purpose: Graph code of the cutoff neighbour list.

Parameters: none

Returns: The integer 0.

Ownership: No allocation.

Errors: None.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `bindings/julia/Seams/test/runtests.jl`

### Seams.MUTUAL

Purpose: Graph code of the mutual four-nearest list.

Parameters: none

Returns: The integer 1.

Ownership: No allocation.

Errors: None.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `bindings/julia/Seams/test/runtests.jl`

### Seams.OVITO

Purpose: CHILL+ rule code used by OVITO.

Parameters: none

Returns: The integer 1.

Ownership: No allocation.

Errors: None.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `bindings/julia/Seams/test/runtests.jl`

### Seams.PUBLISHED

Purpose: CHILL+ rule code of the 2015 paper.

Parameters: none

Returns: The integer 0.

Ownership: No allocation.

Errors: None.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `bindings/julia/Seams/test/runtests.jl`

### Seams.System

Purpose: Empty system from seams_system_create. Positions are added with set_ortho!.

Parameters: none

Returns: A System.

Ownership: The finalizer calls seams_system_destroy.

Errors: Error when seams_system_create returns NULL.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `bindings/julia/Seams/test/runtests.jl`

### Seams.chill_plus

Purpose: CHILL+ labels and the nine-bin census.

Parameters: sys; cutoff; graph; rule

Returns: A pair of vectors, labels and counts.

Ownership: The caller owns the vectors. The System still owns the frame.

Errors: Error with seams_last_error text when the status is not SEAMS_OK.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `bindings/julia/Seams/test/runtests.jl`

### Seams.num_threads

Purpose: Worker threads the engine is using.

Parameters: none

Returns: An Int.

Ownership: No allocation.

Errors: None.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `bindings/julia/Seams/test/runtests.jl`

### Seams.q4q6

Purpose: Local q4 and q6.

Parameters: sys; cutoff

Returns: A pair of vectors, q4 and q6.

Ownership: The caller owns the vectors.

Errors: Error with seams_last_error text when the status is not SEAMS_OK.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `bindings/julia/Seams/test/runtests.jl`

### Seams.set_num_threads

Purpose: Request n worker threads.

Parameters: n

Returns: Nothing.

Ownership: No allocation.

Errors: Error with seams_last_error text when the status is not SEAMS_OK.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `bindings/julia/Seams/test/runtests.jl`

### Seams.set_ortho!

Purpose: Store an orthogonal cell. xyz is 3 by n.

Parameters: sys; xyz; lengths

Returns: The same System.

Ownership: The System still owns the C frame.

Errors: Error when xyz is not 3 by n, or when the C status is not SEAMS_OK.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `bindings/julia/Seams/test/runtests.jl`

### Seams.six_rings

Purpose: Number of primitive six-rings and their fraction.

Parameters: sys; cutoff

Returns: A pair, the six-ring count and the fraction.

Ownership: The temporary ring set is freed with seams_rings_destroy before the function returns.

Errors: Error with seams_last_error text when the status is not SEAMS_OK.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `bindings/julia/Seams/test/runtests.jl`

### Seams.version

Purpose: Engine version string from seams_version_string.

Parameters: none

Returns: A String.

Ownership: Julia owns the string. The C pointer is not freed.

Errors: None.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `bindings/julia/Seams/test/runtests.jl`

### cages

Purpose: Hexagonal and double-diamond cage counts.

Parameters: frame; cutoff

Returns: The Julia value named in the purpose.

Ownership: Returned arrays are owned by Julia. A Frame is freed by its finalizer.

Errors: Error() with the text of seams_last_error when the status is not zero.

Units: cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_julia.jl`

### chill_plus

Purpose: CHILL+ label of each atom.

Parameters: frame; cutoff

Returns: The Julia value named in the purpose.

Ownership: Returned arrays are owned by Julia. A Frame is freed by its finalizer.

Errors: Error() with the text of seams_last_error when the status is not zero.

Units: cutoff in angstrom; labels are dimensionless codes

Example: `docs/examples/minimal_julia.jl`

### ice_ih_path

Purpose: Path of the committed Ice Ih dump, input/traj/ice1h.lammpstrj.

Parameters: none

Returns: The absolute path string.

Ownership: Julia owns the string.

Errors: None.

Units: cutoff and coordinates in angstrom when the call takes them

Example: `docs/examples/minimal_julia.jl`

### jump_tau90

Purpose: Time for an H-H vector to rotate 90 degrees between two frames.

Parameters: first; second; dt

Returns: The Julia value named in the purpose.

Ownership: Returned arrays are owned by Julia. A Frame is freed by its finalizer.

Errors: Error() with the text of seams_last_error when the status is not zero.

Units: dt and the result in the caller's time unit

Example: `docs/examples/minimal_julia.jl`

### neighbour_pairs

Purpose: Directed pairs in the neighbour list.

Parameters: frame; cutoff; type; graph

Returns: The Julia value named in the purpose.

Ownership: Returned arrays are owned by Julia. A Frame is freed by its finalizer.

Errors: Error() with the text of seams_last_error when the status is not zero.

Units: cutoff in angstrom; the count is dimensionless

Example: `docs/examples/minimal_julia.jl`

### nframes

Purpose: Number of frames in a trajectory path.

Parameters: path

Returns: The Julia value named in the purpose.

Ownership: Returned arrays are owned by Julia. A Frame is freed by its finalizer.

Errors: Error() with the text of seams_last_error when the status is not zero.

Units: dimensionless count

Example: `docs/examples/minimal_julia.jl`

### nop

Purpose: Number of atoms in a frame.

Parameters: frame

Returns: The Julia value named in the purpose.

Ownership: Returned arrays are owned by Julia. A Frame is freed by its finalizer.

Errors: Error() with the text of seams_last_error when the status is not zero.

Units: dimensionless count

Example: `docs/examples/minimal_julia.jl`

### primitive_rings

Purpose: Histogram of primitive rings.

Parameters: frame; cutoff; max_ring

Returns: The Julia value named in the purpose.

Ownership: Returned arrays are owned by Julia. A Frame is freed by its finalizer.

Errors: Error() with the text of seams_last_error when the status is not zero.

Units: cutoff in angstrom; bin k is a dimensionless ring count

Example: `docs/examples/minimal_julia.jl`

### q4_q6

Purpose: Local q4 and q6.

Parameters: frame; cutoff

Returns: The Julia value named in the purpose.

Ownership: Returned arrays are owned by Julia. A Frame is freed by its finalizer.

Errors: Error() with the text of seams_last_error when the status is not zero.

Units: cutoff in angstrom; q4 and q6 are dimensionless

Example: `docs/examples/minimal_julia.jl`

### read_frame

Purpose: Read one frame from a dump.

Parameters: path; frame; type

Returns: The Julia value named in the purpose.

Ownership: Returned arrays are owned by Julia. A Frame is freed by its finalizer.

Errors: Error() with the text of seams_last_error when the status is not zero.

Units: frame index is dimensionless; coordinates are angstrom

Example: `docs/examples/minimal_julia.jl`

### signature_cages

Purpose: Closed cages matching a ring-census spec.

Parameters: frame; cutoff; spec

Returns: The Julia value named in the purpose.

Ownership: Returned arrays are owned by Julia. A Frame is freed by its finalizer.

Errors: Error() with the text of seams_last_error when the status is not zero.

Units: cutoff in angstrom; spec is text

Example: `docs/examples/minimal_julia.jl`

### steinhardt

Purpose: Mean local Steinhardt q_l.

Parameters: frame; cutoff; order_l

Returns: The Julia value named in the purpose.

Ownership: Returned arrays are owned by Julia. A Frame is freed by its finalizer.

Errors: Error() with the text of seams_last_error when the status is not zero.

Units: cutoff in angstrom; ql is dimensionless

Example: `docs/examples/minimal_julia.jl`

### symmetry

Purpose: Shell label for Ih, Ic, or stacking disorder.

Parameters: frame

Returns: The Julia value named in the purpose.

Ownership: Returned arrays are owned by Julia. A Frame is freed by its finalizer.

Errors: Error() with the text of seams_last_error when the status is not zero.

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


## Structure still

```{figure} ../tutorials/molecules/julia/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug julia` writes the 2400 by 1800 PNG.
```

