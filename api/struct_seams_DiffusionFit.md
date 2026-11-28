<!-- https://docs.dseams.info/api/struct_seams_DiffusionFit.md -->

# struct seams::DiffusionFit {#struct-seams-diffusionfit}

## Overview {#overview}

Linear Einstein diffusion fit coefficients. [More…](#details-structseams-1-1-diffusion-fit)

## Detailed Documentation {#detailed-documentation}

Linear Einstein diffusion fit coefficients.

### Fields {#fields}

3D diffusion coefficient (slope / 6).

1D axial diffusion coefficient (slope / 2).

2D lateral diffusion coefficient (slope / 4).

Log-log anomalous diffusion exponent, NaN when the fit was skipped.

False when the window had < 2 usable log-log lags.
