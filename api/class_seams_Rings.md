<!-- https://docs.dseams.info/api/class_seams_Rings.md -->

# class seams::Rings {#class-seams-rings}

## Overview {#overview}

Owned primitive-ring set. [More…](#details-classseams-1-1-rings)

### Inherited Members {#inherited-members}

## Detailed Documentation {#detailed-documentation}

Owned primitive-ring set.

Move-only. `begin` and `end` walk the rings and yield vertex indices.

### Methods {#methods}

Number of primitive rings.

**Returns:**

The ring count. A moved-from handle reports 0.

Atom indices of one ring.

**Parameters:**

| index | Ring index, from 0 to size() - 1. |
| --- | --- |

**Returns:**

Frame indices, in ring order.

Histogram of ring size.

**Parameters:**

| up_to | Highest size stored. counts[k] is the number of rings of size k. |
| --- | --- |

**Returns:**

The histogram. Length is `up_to + 1`.

Fraction of atoms that sit on a ring of one size.

**Parameters:**

| ring_size | Size to count. Zero counts every ring. |
| --- | --- |

**Returns:**

Dimensionless fraction.

First ring.

Past the last ring.
