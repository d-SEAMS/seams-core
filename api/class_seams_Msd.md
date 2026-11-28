<!-- https://docs.dseams.info/api/class_seams_Msd.md -->

# class seams::Msd {#class-seams-msd}

## Overview {#overview}

Owned streaming MSD accumulator. [More…](#details-classseams-1-1-msd)

### Inherited Members {#inherited-members}

## Detailed Documentation {#detailed-documentation}

Owned streaming MSD accumulator.

### Construction {#construction}

Allocate a streaming MSD accumulator.

**Parameters:**

| natoms | Particle count per frame. |
| --- | --- |
| max_lag | Maximum lag in frames, or 0 for all lags. |
| origin_stride | Frame stride between origins. |
| dt | Frame spacing in physical time units. |
| axis | Axial direction 3-vector. |

### Methods {#methods}

Push one unwrapped coordinate frame.

Push one periodic coordinate frame and unwrap incrementally.

Copy out accumulated total, axial, and lateral MSD curves.

Fit Einstein diffusion coefficients over `[t_start, t_end]`.
