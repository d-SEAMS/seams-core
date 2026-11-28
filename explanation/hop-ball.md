<!-- https://docs.dseams.info/explanation/hop-ball.md -->

# Carried rings and the hop ball {#carried-rings-and-the-hop-ball}

`primitive::RingUpdater` returns the primitive rings of one neighbour
graph. The list equals `primitive::ringNetwork` on that frame, in the
same canonical order, so the two compare equal element by element.
The measured speedups are the how-to
[Linkcell neighbours and carried updates](../user/howto/linkcell-carried.md).
The primitive test is
[Primitive rings and the full graph](primitive-rings.md).

## What can change {#what-can-change}

Let $G$ and $G'$ be the neighbour graphs of two consecutive frames
with the same vertex count, and let $C$ be the vertices whose
adjacency row differs. Those vertices are the endpoints of every added
or removed edge. Write $D$ for the deepest ring the search keeps, and

$$
L = \lfloor D/2 \rfloor, \qquad
r = \max(0,\ \lfloor (L - 2)/2 \rfloor).
$$

A cycle of at most $D$ members is a primitive ring of one graph and
not of the other only if it contains a changed edge, or a changed edge
lies on a shortcut between two of its members. The shortcut is a path
of $l < d$ edges between two members that are $d$ apart around the
ring, with $d \le L$. One endpoint of that changed edge then lies
within $r$ edges of a ring member, measured along the shortcut.
Every ring that differs therefore meets $T$, the vertices within
distance $r$ of $C$ in the union of $G$ and $G'$.

For $D \le 7$, $r = 0$, so $T = C$. The update drops the stored
rings that meet $T$ and adds the rings of $G'$ that meet $T$. Each
new ring is enumerated once, from its lowest-indexed member in $T$.
When $T$ is every vertex, the search starts at every source.

The first frame, and any frame whose vertex count differs, searches
every source. A frame with no changed row returns the stored list.
For $D \ge 8$ a carried update keeps a ball of radius $L - 1$ at
each vertex and refreshes a ball when an edge within that radius of
the vertex changes. A full search at that depth fills every ball.
For $D \le 7$ the search does not keep those balls.

## When a frame walk carries {#when-a-frame-walk-carries}

A `seams` range keeps one updater when the frames arrive in order on
one worker, `--rebuild` is absent, and `SEAMS_CARRY` is left on. On a
LAMMPS dump the worker count is `--jobs 1`, and the dump is not split
across MPI ranks. On a `.xyz`, `.con`, `.pdb`, `.gro`, or `.dcd` path
the count may be any `--jobs`, because those suffixes are walked by
one worker. A single frame (`--last` omitted or equal to `--frame`),
`--rebuild`, `SEAMS_CARRY=0`, an MPI split of a LAMMPS dump, and a
LAMMPS dump read with `--jobs` above 1 call `ringNetwork` on every
frame. The rings match the carried walk. The flags are in the
[seams CLI](../reference/cli.md) reference.

## How the changed sources run {#how-the-changed-sources-run}

The search over $T$ uses `sourceTeam`. Inside an OpenMP region the
width is one. Outside it, the width stays one while the source count
is below 80 times the OpenMP width. Past that threshold the width is
the source count divided by 80, and at most the OpenMP width.
`teamFor` returns the OpenMP width for any count it is passed.

When `sourceTeam` is one, the new rings are copied and sorted on the
calling thread. When it is wider, the ring indices are sorted with a
stable parallel sort. Outside an OpenMP region the merge of the rings
that miss $T$ uses `teamFor`. Inside a region the merge stays on the
calling thread. Width one walks the stored list and inserts the new
rings in order. A wider width marks the rings that meet $T$ and
merges the rings that stay.

`seams_frame_chill_rings` can keep a cutoff neighbour list across
calls. A positive skin on the cutoff graph builds that list at cutoff
plus skin, and the next call reuses it while the box, cutoff, skin,
type, and graph are unchanged and every atom’s minimum-image
displacement is at most half the skin. The rings of that call are
`ringNetwork`. The ring carry is the updater above, on a sequential
`seams` frame walk.

The ring calls are on [C](../stable/c.md) and [C++](../stable/cpp.md).
