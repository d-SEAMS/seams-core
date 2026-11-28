<!-- https://docs.dseams.info/api/namespace_nucleation.md -->

# namespace nucleation {#namespace-nucleation}

## Overview {#overview}

## Detailed Documentation {#detailed-documentation}

### Global Variables {#global-variables}

Published CHILL+ cutoff. Not a fitted parameter.

Primitive rings through this size. Depth 7 is the adjacency-row search.

Lechner-Dellago $\bar{q}_6$ above this value is an ice-like environment. Set from the gap between a perfect mW Ic crystal and a liquid reference, not from the nucleation trajectory.

q6 solid criterion of ten Wolde, Ruiz-Montero and Frenkel, J. Chem. Phys. 104, 9932 (1996), in the mW form of Haji-Akbari and Debenedetti, PNAS 112(34) (2015), doi:10.1073/pnas.1509267112, SI Table S1: neighbours within 3.2 A (the first minimum of the mW g(r)), and a molecule is solid when the mean over its neighbours of the normalised q6(i).q6*(j) exceeds 0.5. The solid cluster is the largest component of solid molecules on the same 3.2 A graph. The Reinhardt chain exclusion is not applied.

### Global Functions {#global-functions}

Per-molecule ten Wolde q6 score: the mean over neighbours within cutoff of `Re q6(i).q6*(j) / (|q6(i)| |q6(j)|)`. A molecule without neighbours scores 0. The box is orthorhombic.

Label one frame. Leaves iceType as the published CHILL+ label. typeI selects the molecules passed to q6-bar. q6Cut is the ice cut on $\bar{q}_6$.

Cubic / (cubic + hexagonal). Quiet NaN when the denominator is 0.

Ring signal is n6 / (n5 + n6). Ice signal is the ice-like count divided by the molecule count. Both are intensive.

How far definition b moves the onset size and the onset frame relative to definition a. dSize = b.size - a.size at each onset.
