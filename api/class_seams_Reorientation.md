<!-- https://docs.dseams.info/api/class_seams_Reorientation.md -->

# class seams::Reorientation {#class-seams-reorientation}

## Overview {#overview}

Owned streaming orientational correlation accumulator. [More…](#details-classseams-1-1-reorientation)

### Inherited Members {#inherited-members}

## Detailed Documentation {#detailed-documentation}

Owned streaming orientational correlation accumulator.

### Construction {#construction}

Allocate a reorientation accumulator.

**Parameters:**

| n_vectors | Molecular vector count per frame. |
| --- | --- |
| max_lag | Maximum lag in frames, or 0 for all lags. |
| origin_stride | Frame stride between origins. |
| dt | Frame spacing. |

### Methods {#methods}

Push one frame of molecular 3-vectors.

Read C1 and C2 correlation curves.

Fit integrated and single-exponential reorientation times.
