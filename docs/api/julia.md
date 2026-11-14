# Julia

`bindings/julia/Dseams` is the package CI runs. It calls `libseams_lang`. Point `SEAMS_LANG_LIB` at that shared library and `DSEAMS_ROOT` at the checkout. `bindings/julia/Seams` is the thin `libseams` example.

`meson test docs_api_julia` runs `docs/examples/julia_api.jl` with
`SEAMS_LANG_LIB` pointed at `libseams_lang`.
`docs/examples/minimal_julia.jl` reads the
GenIce sI cell and calls `steinhardt`, `chill_plus`, and `primitive_rings`.
`meson test docs_examples` runs it.

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
<!-- API_SYMBOLS_END -->


## Structure still

```{figure} ../tutorials/molecules/julia/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug julia` writes the 2400 by 1800 PNG.
```

