# Primitive and chordless rings

Lean 4 development of the ring facts used by the d-SEAMS search: primitive
rings are chordless, the converse fails on an explicit graph, enumeration
from shortest-path pairs is complete, pairs on those paths need not be
tested again, a shortcut on a cycle of length at most 7 is a walk of length
at most 2, and a bond change moves primitive rings only inside a ball of
radius `k / 2`.

Build with `lake exe cache get` followed by `lake build`.
