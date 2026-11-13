# Science front

These calls live in `src/include/internal/science_front.hpp` on
`science/interface-nucleation` (`55192e3d`). They are not in `seams.h`.
A build of this docs branch does not compile them. `science-front` is
the driver in `scripts/science_front.cpp`.

The cutoff is `science::kCutoffA`, 3.5 Å. `science::kQ6BarIceCut` is
0.36. Neither value is fitted. CHILL+ is the published Nguyen and
Molinero rule: q3 over the four nearest neighbours, a crystal row of
four neighbours inside the cutoff, and the interfacial neighbour test.
The OVITO rule is a separate label. It does not choose these labels.
The shell label is 0 hexagonal, 1 cubic, 2 stacking-disordered, and
-1 otherwise. q6-bar is the Lechner-Dellago average on the same cutoff.

## Layers

`science::layersOf(cloud, axis, binWidth)` bins CHILL+, the shell
label, q6 and q6-bar along one box edge. `axis` is 0, 1 or 2. A
negative axis uses the longest edge. The default bin width is 1 Å.
Empty bins between ice planes are filled before the front is fit.
`Layers::n` is the number of molecules. `Layers::rows` is one
`LayerRow` per bin: the bin centre, the occupancy, the cubic,
hexagonal, interfacial and liquid counts, the shell counts, and the
mean q6 and q6-bar.

`science::fitSolidFraction(phi, mid, span)` turns a binned solid
fraction into the 10-90 fronts. `Front::position` is the 0.5 crossing.
`Front::width` is the 10-90 width of that crossing. A flat crystal,
with no fall through 0.5, leaves `Layers::defined` false and the front
list empty.

`science::bootstrapSlope(x, y, nDraw, seed)` is the ordinary-least-squares
slope and a percentile bootstrap. A draw whose x values have no spread
is dropped. The driver uses it for the two front positions and for the
ice thickness. Half the thickness slope is the mean advance of the two
fronts when the slab grows from both sides.

## One frame

`science::nucleusOf(cloud, frame, q6Cut)` labels one frame and leaves
the published CHILL+ label on the cloud. `NucleusFrame` holds the
CHILL+ census, the shell census, the q6-bar ice count and mean, the
stacking molecules split by the CHILL+ label, and the polytype
agreement. `core`, `extended`, `shell` and `q6` are the largest
clusters. The extended cluster is cubic, hexagonal and interfacial.
The core cluster is cubic and hexagonal.

`science::labelAtoms(cloud, q6Cut)` returns one `AtomRow` per molecule:
coordinates, the published CHILL+ label, the OVITO CHILL+ label, the
shell label, q6, q6-bar, and whether the molecule is in the extended
cluster.

## Critical size

`science::criticalSize(series, nDraw, seed)` reads the largest-cluster
size at each frame. The critical size is the mean of the contiguous
band around the smallest value. The band stops where the size exceeds
that minimum by the larger of 5 molecules and 15 percent. `nStar` is
that mean. `lo` and `hi` are a percentile bootstrap of the mean.
`onset` is the frame of the minimum. `window` is the width of the band.

This is an upper bound, not a committor. A minimum in the first two
frames, followed by a larger cluster, sets `upperBound`. The dump
starts after the nucleus is already at least that large. The number
is not the probability that the cluster will grow.
