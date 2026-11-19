<!-- https://docs.dseams.info/api/class_seams_Fingerprint.md -->

# class seams::Fingerprint {#class-seams-fingerprint}

## Overview {#overview}

Owned topological fingerprint. [More…](#details-classseams-1-1-fingerprint)

### Inherited Members {#inherited-members}

## Detailed Documentation {#detailed-documentation}

Owned topological fingerprint.

### Methods {#methods}

Frame fingerprint text.

**Returns:**

The key, without a trailing NUL.

Local key of one atom.

**Parameters:**

| index | Atom index in the frame. |
| --- | --- |

**Returns:**

The local key.

Primitive-ring census stored on the fingerprint.

**Returns:**

`counts[k]` is the number of rings of size k. Length 13.
