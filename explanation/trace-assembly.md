<!-- https://docs.dseams.info/explanation/trace-assembly.md -->

# TRACE and signature cages {#trace-and-signature-cages}

TRACE (Hsu and Lin, J. Chem. Theory Comput. 2025,
[10.1021/acs.jctc.5c01459](https://doi.org/10.1021/acs.jctc.5c01459))
finds hydrate cages on a hydrogen-bond graph. d-SEAMS answers the oxygen
question: which primitive rings tile a closed polyhedron of a given face
census. The timed comparison is the how-to
[Cages and TRACE](../user/howto/cages-trace.md).

## The hydrogen-bond cage {#the-hydrogen-bond-cage}

A hydrogen bond is kept when the donor-acceptor oxygen distance is below
a cutoff and the angle at the donor oxygen, between the donated hydrogen
and the acceptor oxygen, is within the angle cutoff of a linear bond.
Rings are enumerated by depth-first search. A cup is a base ring plus
the lateral rings that share one edge with it. The cup closes when the
last lateral ring meets both of its neighbours. Complete and incomplete
cages are classified from those cups.

The TRACE program’s default cutoffs are 0.36 nm and 35°. The comparison
calls ClavisCage 0.3.18 at 0.35 nm and 35°. On an oxygen-only frame that
comparison places one hydrogen 0.96 Å from each oxygen along each of two
ice-rule bonds, then hands the waters to TRACE. d-SEAMS does not place
those hydrogens and does not call TRACE.

The three censuses in the comparison are $5^{12}$, $5^{12}6^{2}$ and
$5^{12}6^{4}$. On ice, which has no five-rings, both searches report
no cage of those three kinds. The counts are in the how-to above.

## A closed face list {#a-closed-face-list}

`seams_frame_signature` and `seams cages --signature` count a set of
primitive rings as one cage when three conditions hold. Every edge of
every face is shared by exactly two faces of the set. At least three
faces around each vertex form one cycle. The cell is a topological
sphere:

$$
V - E + F = 2.
$$

$F$ is the number of faces, $E$ the number of distinct edges, and
$V$ the number of distinct vertices. A ring whose minimum-image bond
shifts sum to a nonzero lattice vector winds through the periodic cell
and is dropped before this test. A cell thinner than twice the neighbour
cutoff is refused, because one pair can then meet through more than one
periodic image.

The default face-edge window is $[2.65, 3.15]$ in the coordinate units
of the file. `none` leaves the bound to the face census alone. On the
46-molecule sI cell, `--signature 512` reports 2 cages and 40 atoms with
the default window and with `none`.

The string `512` is twelve pentagons. `51262` is twelve pentagons and
two hexagons. `51264` is twelve pentagons and four hexagons. Joining
them with `+` prints one line per census. The call that takes the window
is `seams_frame_signature_edge`. The C API names are
`seams_signature_cages` and `seams_signature_cages_edge`.

[C++](../stable/cpp.md) documents `seams_frame_signature` and
`seams_frame_signature_edge`.
[C](../stable/c.md) documents `seams_signature_cages` and
`seams_signature_cages_edge`.
