# srr standards assessment

`dseams` is submitted to rOpenSci software peer review. It is not
statistical software. The [srr standards](https://stats-devguide.ropensci.org/standards.html)
(version 0.2.0) apply to packages in a statistical category. This note
records which of those standards are relevant and why the package does not
carry `@srrstats` tags.

## Categories that do not apply

The category standards are for software whose purpose is statistical
inference or statistical modeling:

- Bayesian
- regression and supervised learning
- machine learning
- probability distributions
- exploratory data analysis
- spatial
- time series
- unsupervised learning

`dseams` labels local structure in a molecular configuration (CHILL+,
primitive rings, cages, Steinhardt parameters, shell symmetry). It does not
fit or sample a statistical model, and it does not return standard errors,
posteriors, or predictions from a trained model. Those category standards
are not applicable. Adding `@srrstats` tags would route the package onto
the statistical-software checklist, which is the wrong review track.

## General standards that overlap this package

The general standards are written for statistical software. A few of them
describe ordinary scientific-software practice, and the package already
does that work:

| Standard | Where it is met |
| --- | --- |
| G1.0 Primary published reference | DESCRIPTION and the roxygen `@references` tags cite Goswami, Goswami and Singh (2020), Nguyen and Molinero (2015), Franzblau (1991), and Steinhardt, Nelson and Ronchetti (1983) |
| G1.1 What the algorithm is | The functions implement those published deterministic classifiers. They are not a new statistical estimator |
| G1.2 Life cycle | Version 2.14.0. `NEWS.md` records this as the first CRAN submission of a stable engine binding |
| G1.4 roxygen2 | Every exported function is documented with roxygen2, and each has an example |
| G5.0, G5.1, G5.4 Known data | `inst/extdata` ships the ice Ih cell. Tests require 128 hexagonal waters and 256 primitive six-rings |
| G5.2 Error behavior | Tests cover a missing file, a file that is not a dump, an unknown CHILL+ rule, an unknown example name, and a non-positive thread count |
| G5.3 Finite results | The ice Ih test requires finite `q6` |
| G2.0, G2.1 Input checks | `dseams_threads()` rejects a non-positive count. `dseams_example()` rejects an unknown name. `jump_tau90()` requires `dseams_frame` objects |

## General standards that are not relevant

These assume a statistical implementation. They are not applicable, and
they are not tagged `@srrstatsNA` for the reason above.

- G1.3, G1.5, G1.6: statistical terminology, and performance claims against other statistical packages. The package does not make statistical performance claims.
- G2.2 through G2.16: tabular statistical input, factors, and missing-data policies for model frames. Inputs are trajectory paths and a `dseams_frame`, not a model matrix.
- G3.0 through G3.1a: floating-point equality in statistical algorithms, and choices of covariance estimator. The R code does not compare floating-point results for equality. Covariance is not computed.
- G4.0: file-suffix handling for statistical output files. The package does not write result files.
- G5.5 through G5.12: random seeds, parameter recovery, noise susceptibility, and extended statistical test batteries. The classifiers are deterministic. The reference result is the published ice Ih cell, not a simulation draw.
