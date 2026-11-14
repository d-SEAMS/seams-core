# dseams

Ice classification and ring statistics for molecular simulations. `dseams`
is the R interface to the [d-SEAMS](https://dseams.info) engine: CHILL+
labels, primitive rings, hexagonal and double-diamond cages, and Steinhardt
order parameters.

## Installation

The development version is the `bindings/r/dseams` directory of
[seams-core](https://github.com/HaoZeke/seams-core). `configure` copies the
engine sources from that checkout, so the install does not download them.

```r
# install.packages("remotes")
remotes::install_github("HaoZeke/seams-core", subdir = "bindings/r/dseams")
```

From a seams-core checkout:

```sh
R CMD INSTALL bindings/r/dseams
```

A C++20 compiler is required (`SystemRequirements` in `DESCRIPTION`).

## Use

```r
library(dseams)
frame <- read_frame(dseams_example("ice1h.lammpstrj"), type = -1L)
table(chill_plus(frame)$name)
primitive_rings(frame, graph = "cutoff", max_depth = 6L)
```

`vignette("classifying-ice", package = "dseams")` walks through the same
dumps. Report bugs at <https://github.com/d-SEAMS/seams-core/issues>.

The package is MIT licensed. Cite Goswami, Goswami and Singh (2020)
<doi:10.1021/acs.jcim.0c00031>. `citation("dseams")` prints the record.
