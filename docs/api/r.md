# R

The package is `bindings/r/dseams`. `libseams_lang` serves the frame
generics. `libseams` serves `seams_system_new` and the analyses named
`seams_*`. roxygen2 writes `man/` and `NAMESPACE`. pkgdown renders the
package site from `_pkgdown.yml`.

```bash
meson setup build -Dwith_tests=true
meson compile -C build
DSEAMS_ROOT=$PWD R CMD INSTALL bindings/r/dseams
```

`R CMD check` is the cleanliness check. It expects no errors and no
warnings. `pkgdown::build_site()` writes the native HTML.

## Quickstart

```r
library(dseams)
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
:lines: 9-36,77-91
```

<!-- API_SYMBOLS_BEGIN -->

## Symbol reference

### cages

Purpose: Cage counts for a dump path.

Parameters: path; frame; cutoff

Units: cutoff in angstrom; frame index is dimensionless

Example: `docs/examples/minimal_r.R`

### cages.dseams_frame

Purpose: Cage counts for a frame object.

Parameters: frame; cutoff

Units: cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_r.R`

### chill_plus

Purpose: CHILL+ labels for a dump path.

Parameters: path; frame; cutoff

Units: cutoff in angstrom

Example: `docs/examples/minimal_r.R`

### chill_plus.dseams_frame

Purpose: CHILL+ labels for a frame object.

Parameters: frame; cutoff

Units: cutoff in angstrom; labels are dimensionless codes

Example: `docs/examples/minimal_r.R`

### jump_tau90

Purpose: Rotation time of H-H vectors between two frames.

Parameters: first; second; dt

Units: dt and the result in the caller's time unit

Example: `docs/examples/minimal_r.R`

### nframes

Purpose: Number of frames in a trajectory path.

Parameters: path

Units: dimensionless count

Example: `docs/examples/minimal_r.R`

### positions

Purpose: Coordinates of a dump frame.

Parameters: path; frame

Units: angstrom

Example: `docs/examples/minimal_r.R`

### positions.dseams_frame

Purpose: Coordinates stored on a frame object.

Parameters: frame

Units: angstrom

Example: `docs/examples/minimal_r.R`

### primitive_rings

Purpose: Ring histogram for a dump path.

Parameters: path; frame; cutoff

Units: cutoff in angstrom

Example: `docs/examples/minimal_r.R`

### primitive_rings.dseams_frame

Purpose: Ring histogram for a frame object.

Parameters: frame; cutoff

Units: cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_r.R`

### print.dseams_frame

Purpose: Print the atom count of a frame.

Parameters: frame

Units: dimensionless count

Example: `docs/examples/minimal_r.R`

### q4_q6

Purpose: Local q4 and q6 for a dump path.

Parameters: path; frame; cutoff

Units: cutoff in angstrom; q4 and q6 are dimensionless

Example: `docs/examples/minimal_r.R`

### q4_q6.dseams_frame

Purpose: Local q4 and q6 for a frame object.

Parameters: frame; cutoff

Units: cutoff in angstrom; q4 and q6 are dimensionless

Example: `docs/examples/minimal_r.R`

### read_frame

Purpose: Read one frame from a dump.

Parameters: path; frame; type

Units: frame index is dimensionless; coordinates are angstrom

Example: `docs/examples/minimal_r.R`

### signature_cages

Purpose: Closed cages for a dump path.

Parameters: path; frame; spec

Units: spec is text; cutoff used by the call is in angstrom

Example: `docs/examples/minimal_r.R`

### signature_cages.dseams_frame

Purpose: Closed cages for a frame object.

Parameters: frame; spec

Units: spec is text

Example: `docs/examples/minimal_r.R`

### steinhardt

Purpose: Mean local q_l for a dump path.

Parameters: path; frame; order_l

Units: order parameter is dimensionless

Example: `docs/examples/minimal_r.R`

### steinhardt.dseams_frame

Purpose: Mean local q_l for a frame object.

Parameters: frame; order_l

Units: order parameter is dimensionless

Example: `docs/examples/minimal_r.R`

### symmetry

Purpose: Shell labels for a dump path.

Parameters: path; frame

Units: dimensionless label codes

Example: `docs/examples/minimal_r.R`

### symmetry.dseams_frame

Purpose: Shell labels for a frame object.

Parameters: frame

Units: dimensionless label codes

Example: `docs/examples/minimal_r.R`
### seams_system_new

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_system_read

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_system_count

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_chill

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_steinhardt

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_cages

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_signature_cages

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_neighbours_knn

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_rdf

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_running_cn

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_coordination

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_hbonds

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_density_axis

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_fingerprint

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_contact_pairs

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_domains

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_prisms

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_monolayer

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_rings_vertices

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_site_table

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

<!-- API_SYMBOLS_END -->
