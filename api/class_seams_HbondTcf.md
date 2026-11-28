<!-- https://docs.dseams.info/api/class_seams_HbondTcf.md -->

# class seams::HbondTcf {#class-seams-hbondtcf}

## Overview {#overview}

Owned streaming hydrogen-bond lifetime accumulator. [More…](#details-classseams-1-1-hbond-tcf)

### Inherited Members {#inherited-members}

## Detailed Documentation {#detailed-documentation}

Owned streaming hydrogen-bond lifetime accumulator.

### Construction {#construction}

Allocate an H-bond lifetime accumulator.

**Parameters:**

| max_lag | Maximum lag in frames, or 0 for all lags. |
| --- | --- |
| origin_stride | Frame stride between origins. |
| dt | Frame spacing. |
| directed | True for directed (donor, acceptor) edges. |

### Methods {#methods}

Push one frame’s `(donor, acceptor)` arrays.

Read intermittent C_HB and continuous S_HB curves.

Fit integrated and single-exponential lifetimes.
