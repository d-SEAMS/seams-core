# R

d-SEAMS 2.14.0. The package is `bindings/r/dseams`. `libseams_lang`
serves the frame generics through cpp11. `libseams` serves
`seams_system_new` and the analyses named `seams_*`. roxygen2 writes
`man/` and `NAMESPACE`. pkgdown renders the package site from
`_pkgdown.yml`. The exported frames are `read_frame`, `nframes`,
`positions`, `chill_plus`, `cages`, `primitive_rings`,
`signature_cages`, `symmetry`, `steinhardt`, `q4_q6`, and `jump_tau90`.
`print` on a frame dispatches to `print.dseams_frame`.
The methods are `positions.dseams_frame`, `chill_plus.dseams_frame`,
`cages.dseams_frame`, `primitive_rings.dseams_frame`,
`signature_cages.dseams_frame`, `symmetry.dseams_frame`,
`steinhardt.dseams_frame`, and `q4_q6.dseams_frame`.

`docs/examples/core_counts.R` checks Ice Ih (128 hexagonal, 256
six-rings) and sII (136 clathrate, 144 five-rings, 16 six-rings).
`meson test docs_core_counts` runs it.

```{literalinclude} ../../../docs/examples/core_counts.R
:language: r
```

`dseams` is not on CRAN yet and installs from this checkout. With only
`DSEAMS_ROOT` set, `configure` runs `tools/vendor.sh` and the package
compiles its own copy of the engine. To install the package against a
built library instead, also name `libseams_lang` and `seams_lang.h`:

```bash
meson setup build -Dwith_tests=true
meson compile -C build
DSEAMS_ROOT=$PWD R CMD INSTALL bindings/r/dseams
SEAMS_LIBDIR=$PWD/build/bindings SEAMS_INCLUDE=$PWD/bindings/cxx R CMD INSTALL --preclean bindings/r/dseams
```

`--preclean` clears the objects the first install left in `src`, which
`R CMD INSTALL` would otherwise reuse. `R CMD build bindings/r/dseams` makes a source tarball that carries the
engine and installs with no checkout; `R CMD check --as-cran` on it
leaves the one NOTE for a new submission.

CI runs the testthat file when `Rscript` is on `PATH` and
`DSEAMS_ROOT` points at this tree. The meson target is `language_fronts`
(`tests/run_language_fronts.sh`). A missing R skips that front.
`docs/examples/minimal_r.R` reads the GenIce sI
cell and calls `steinhardt`, `chill_plus`, and `primitive_rings`.
`meson test docs_examples` runs it.

`R CMD check` is the cleanliness check. It expects no errors and no
warnings. `pkgdown::build_site()` writes the native HTML from the
roxygen topics. The built site is linked from this book:

```{raw} html
<p><a href="../r/index.html">dseams pkgdown reference</a>,
generated from the roxygen topics in <code>man/</code>.</p>
```

## Quickstart

```r
library(dseams)
path <- system.file("extdata", "exampleTraj.lammpstrj", package = "dseams")
frame <- read_frame(path, type = 2L)
ice <- chill_plus(frame)
print(ice)
summary(frame)
```

On `input/traj/exampleTraj.lammpstrj`, oxygen type 2, the cutoff shell
has 9 interClathrate and 241 water. The mutual shell has 1 clathrate,
7 interClathrate, and 242 water. Mean q3 is 0.374897. Six-rings at depth
8 number 224.

## S3 results

`read_frame` returns class `dseams_frame`. `chill_plus` returns
`dseams_ice`, `cages` returns `dseams_cages`, `primitive_rings` returns
`dseams_rings`, `steinhardt` returns `dseams_steinhardt`, and `q4_q6`
returns `dseams_q4_q6`. Each has `print` and `summary`. Cage counts are
attributes `n_ih`, `n_ic`, `n_both`, and `n_water`, and they match
`seams cages`.

`seams_system_new` returns an external pointer of class `dseams_system`.
The finalizer calls `seams_system_destroy`. `print` reports the atom
count. Analysis helpers return the same vectors the C oracle checks,
with a class so `print` can name the columns. `as.numeric` and `[[`
keep the positions those tests use.

## Errors

A missing path stops in `normalizePath` before the library is called.
An empty system stops with the text from `seams_handle_error`. The
condition is an ordinary R error, so `tryCatch` and testthat
`expect_error` see it.

## Ownership

The frame pointer inside `dseams_frame` is an external pointer freed by
the C finalizer. Do not keep a copy of `frame$ptr` after the frame is
collected. System, neighbour, ring, cage, and fingerprint handles used
inside one `seams_*` helper are freed before that helper returns, except
the system itself, which the caller created and the finalizer releases.

```{literalinclude} ../../../docs/examples/minimal_r.R
:language: r
```

```{literalinclude} ../../../docs/examples/r_api.R
:language: r
```

```{literalinclude} ../../../bindings/r/dseams/tests/testthat/test-ice.R
:language: r
:lines: 22-49,84-102
```

<!-- API_SYMBOLS_BEGIN -->

## Symbol reference

### cages

Purpose: Cage counts for a dump path.

Parameters: path; frame; cutoff

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: cutoff in angstrom; frame index is dimensionless

Example: `docs/examples/minimal_r.R`

### cages.dseams_frame

Purpose: Cage counts for a frame object.

Parameters: frame; cutoff

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_r.R`

### chill_plus

Purpose: CHILL+ labels for a dump path.

Parameters: path; frame; cutoff

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: cutoff in angstrom

Example: `docs/examples/minimal_r.R`

### chill_plus.dseams_frame

Purpose: CHILL+ labels for a frame object.

Parameters: frame; cutoff

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: cutoff in angstrom; labels are dimensionless codes

Example: `docs/examples/minimal_r.R`

### dseams_example

Purpose: Path of a trajectory shipped in the package's `inst/extdata`: `ice1h.lammpstrj`, `exampleTraj.lammpstrj` or `nucleation_f1.lammpstrj`. With no name it lists them.

Parameters: name (file name, or `NULL` to list the files)

Returns: A file path, or the file names when name is NULL.

Ownership: R owns the returned character vector.

Errors: The binding stops with "no example named" when name is not one of the shipped files.

Units: none; returns a file path or file names

Example: `docs/examples/minimal_r.R`

### dseams_threads

Purpose: Report the engine thread count, or set it for later calls from this R session. Under `_R_CHECK_LIMIT_CORES_` the package starts with at most two threads; the option `dseams.threads` sets the starting count.

Parameters: n (number of threads, or `NULL` to only query)

Returns: The thread count in effect after the call.

Ownership: No allocation.

Errors: The binding stops when n is not one positive integer.

Units: dimensionless thread count; 1 when the package was built without OpenMP

Example: `docs/examples/minimal_r.R`

### dseams_version

Purpose: Version string of the seams-core engine the package was built with or linked against.

Parameters: none

Returns: The engine version string.

Ownership: R owns the returned character vector.

Errors: None.

Units: none; returns a version string

Example: `docs/examples/minimal_r.R`

### jump_tau90

Purpose: Rotation time of H-H vectors between two frames.

Parameters: first; second; dt

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: dt and the result in the caller's time unit

Example: `docs/examples/minimal_r.R`

### nframes

Purpose: Number of frames in a trajectory path.

Parameters: path

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: dimensionless count

Example: `docs/examples/minimal_r.R`

### positions

Purpose: Coordinates of a dump frame.

Parameters: path; frame

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom

Example: `docs/examples/minimal_r.R`

### positions.dseams_frame

Purpose: Coordinates stored on a frame object.

Parameters: frame

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom

Example: `docs/examples/minimal_r.R`

### primitive_rings

Purpose: Ring histogram for a dump path.

Parameters: path; frame; cutoff

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: cutoff in angstrom

Example: `docs/examples/minimal_r.R`

### primitive_rings.dseams_frame

Purpose: Ring histogram for a frame object.

Parameters: frame; cutoff

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_r.R`

### print.dseams_cages

Purpose: Print cage counts.

Parameters: x; ...

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: as named by the parameters

Example: `docs/examples/minimal_r.R`

### print.dseams_chill

Purpose: Print CHILL population counts.

Parameters: x; ...

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: as named by the parameters

Example: `docs/examples/minimal_r.R`

### print.dseams_frame

Purpose: Print the atom count of a frame.

Parameters: frame

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: dimensionless count

Example: `docs/examples/minimal_r.R`

### print.dseams_ice

Purpose: Print CHILL+ labels.

Parameters: x; ...

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: as named by the parameters

Example: `docs/examples/minimal_r.R`

### print.dseams_q4_q6

Purpose: Print q4 and q6.

Parameters: x; ...

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: as named by the parameters

Example: `docs/examples/minimal_r.R`

### print.dseams_rings

Purpose: Print a ring histogram.

Parameters: x; ...

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: as named by the parameters

Example: `docs/examples/minimal_r.R`

### print.dseams_steinhardt

Purpose: Print Steinhardt parameters.

Parameters: x; ...

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: as named by the parameters

Example: `docs/examples/minimal_r.R`

### print.dseams_system

Purpose: Print a system handle.

Parameters: x; ...

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: as named by the parameters

Example: `docs/examples/minimal_r.R`

### q4_q6

Purpose: Local q4 and q6 for a dump path.

Parameters: path; frame; cutoff

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: cutoff in angstrom; q4 and q6 are dimensionless

Example: `docs/examples/minimal_r.R`

### q4_q6.dseams_frame

Purpose: Local q4 and q6 for a frame object.

Parameters: frame; cutoff

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: cutoff in angstrom; q4 and q6 are dimensionless

Example: `docs/examples/minimal_r.R`

### read_frame

Purpose: Read one frame from a dump.

Parameters: path; frame; type

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: frame index is dimensionless; coordinates are angstrom

Example: `docs/examples/minimal_r.R`

### seams_cages

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_chill

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_contact_pairs

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_coordination

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_density_axis

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_domains

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_fingerprint

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_hbonds

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_monolayer

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_neighbours_knn

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_prisms

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_rdf

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_rings_vertices

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_running_cn

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_set_num_threads

Purpose: Set the OpenMP width used by the next analysis.

Parameters: n

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: as named by the parameters

Example: `docs/examples/minimal_r.R`

### seams_signature_cages

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_site_table

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_steinhardt

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_system_count

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_system_new

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_system_read

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### signature_cages

Purpose: Closed cages for a dump path.

Parameters: path; frame; spec

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: spec is text; cutoff used by the call is in angstrom

Example: `docs/examples/minimal_r.R`

### signature_cages.dseams_frame

Purpose: Closed cages for a frame object.

Parameters: frame; spec

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: spec is text

Example: `docs/examples/minimal_r.R`

### steinhardt

Purpose: Mean local q_l for a dump path.

Parameters: path; frame; order_l

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: order parameter is dimensionless

Example: `docs/examples/minimal_r.R`

### steinhardt.dseams_frame

Purpose: Mean local q_l for a frame object.

Parameters: frame; order_l

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: order parameter is dimensionless

Example: `docs/examples/minimal_r.R`

### summary.dseams_cages

Purpose: Summary of cage names.

Parameters: object; ...

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: as named by the parameters

Example: `docs/examples/minimal_r.R`

### summary.dseams_chill

Purpose: Summary of CHILL population counts.

Parameters: object; ...

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: as named by the parameters

Example: `docs/examples/minimal_r.R`

### summary.dseams_frame

Purpose: Summary of one frame.

Parameters: object; ...

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: as named by the parameters

Example: `docs/examples/minimal_r.R`

### summary.dseams_ice

Purpose: Summary of CHILL+ labels.

Parameters: object; ...

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: as named by the parameters

Example: `docs/examples/minimal_r.R`

### summary.dseams_q4_q6

Purpose: Summary of q4 and q6.

Parameters: object; ...

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: as named by the parameters

Example: `docs/examples/minimal_r.R`

### summary.dseams_rings

Purpose: Summary of a ring histogram.

Parameters: object; ...

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: as named by the parameters

Example: `docs/examples/minimal_r.R`

### summary.dseams_steinhardt

Purpose: Summary of Steinhardt parameters.

Parameters: object; ...

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: as named by the parameters

Example: `docs/examples/minimal_r.R`

### summary.dseams_system

Purpose: Summary of a system handle.

Parameters: object; ...

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: as named by the parameters

Example: `docs/examples/minimal_r.R`

### symmetry

Purpose: Shell labels for a dump path.

Parameters: path; frame

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: dimensionless label codes

Example: `docs/examples/minimal_r.R`

### symmetry.dseams_frame

Purpose: Shell labels for a frame object.

Parameters: frame

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: dimensionless label codes

Example: `docs/examples/minimal_r.R`
<!-- API_SYMBOLS_END -->


## Structure still

```{figure} ../tutorials/molecules/r/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug r` writes the 2400 by 1800 PNG.
```

