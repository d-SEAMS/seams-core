# R

d-SEAMS 2.14.0. The R package is `bindings/r/dseams`. It calls
`libseams_lang` through cpp11. The exported frames are `read_frame`,
`nframes`, `positions`, `chill_plus`, `cages`, `primitive_rings`,
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

Build the library, then install the package against it:

```bash
meson setup build -Dwith_tests=true
meson compile -C build
DSEAMS_ROOT=$PWD R CMD INSTALL bindings/r/dseams
```

CI runs the testthat file when `Rscript` is on `PATH` and
`DSEAMS_ROOT` points at this tree. The meson target is `language_fronts`
(`tests/run_language_fronts.sh`). A missing R skips that front.
`docs/examples/minimal_r.R` reads the GenIce sI
cell and calls `steinhardt`, `chill_plus`, and `primitive_rings`.
`meson test docs_examples` runs it.

```{literalinclude} ../../../docs/examples/minimal_r.R
:language: r
```
`meson test docs_api_r` installs the package against `libseams_lang` and
runs the script below.

```{literalinclude} ../../../docs/examples/r_api.R
:language: r
```

`seams_lang_fixtures` is the C++ driver of `seams_lang.h`, not the R package.

```bash
DSEAMS_ROOT=$PWD meson test -C build language_fronts --print-errorlogs
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
<!-- API_SYMBOLS_END -->


## Structure still

```{figure} ../tutorials/molecules/r/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug r` writes the 2400 by 1800 PNG.
```

