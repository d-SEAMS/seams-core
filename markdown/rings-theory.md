# Rings and cages: SQQ against the primitive-ring search

The comparison is Pang and Sun, *J. Phys. Chem. A* **130** (41), 8187–8199 (2026), doi:10.1021/acs.jpca.6c05769. The engine read here is SQQ 0.5.8 (`find_rings` and `find_cages` in `sqq/core/sqq_cpp/native/src/core.cpp`, and the same contract in `docs/design.md`). d-SEAMS enumerates Franzblau primitive rings. SQQ’s default faces are chordless cycles. Those are not the same set.

## Two ring predicates

A ring is a simple cycle. Write `ringDist(i, j)` for the shorter of the two arcs.

- **Chordless.** The only edges of the graph between vertices of the cycle are the cycle edges. A chord is an extra edge, nothing longer.
- **Primitive (Franzblau, shortest-path).** For every pair of vertices on the cycle, the graph distance equals `ringDist`. Every primitive ring is chordless, because a chord is a path of length 1 between vertices at ring distance at least 2. The converse is false.

Counterexample. Vertices `0–1–2–3–4–5–0`, plus a vertex `6` joined only to `0` and `3`. The hexagon has no chord. The graph distance from `0` to `3` is 2 and the ring distance is 3, so the hexagon is not primitive. The two pentagons `0–1–2–3–6` and `0–5–4–3–6` are primitive. SQQ’s default search keeps the hexagon. `primitive::ringNetwork` does not. Both facts are checked in `tests/test_ring_hop.cpp`, and the separation is proved in `formal/Formal/Rings.lean` (`hex_chordless`, `hex_not_primitive`, `primitive_chordless`).

A second, smaller split: the square `0–1–2–3` with the diagonal `0–2`. Neither predicate keeps the square. Both keep the two triangles.

SQQ can be asked for `ring.definition: shortest_path`. That filter is Franzblau’s test applied after the chordless search, so it only deletes rings. It does not add any. On a hydrate network the deleted rings change the faces that GROW is allowed to use.

## What each search does

### SQQ rings

`find_rings` (C++, 0.5.8) walks every simple path whose minimum vertex is the start and whose second vertex is less than its last. At each step a neighbour that touches any earlier vertex except the start is a chord and the branch dies. A neighbour that touches the start closes the cycle and the path is not extended past that edge. Closed paths whose length is in the requested set are stored in a `std::set` after a full rotation-and-reversal canonicalisation.

The Python search is the same walk, with an optional second test: for `shortest_path`, a bounded BFS from each ring vertex out to `floor(size/2)` must reproduce the ring distances. Those distance maps are cached per `(source, depth)`.

### SQQ cages

GROW starts from each canonical seed ring. A branch may not add a ring of lower rank than the seed, so a shell is not rediscovered from every face. The state is the sorted face list, the edges used once, and a mask of Euler compositions still possible. An edge used twice is closed. A third use is rejected. The open edge with the fewest eligible rings is expanded, and every eligible ring is tried. Closure requires every edge twice, `V − E + F = 2`, a connected face-adjacency graph, one cyclic link at each vertex, trivalence, and a face count in the generated family `2 n4 + n5 = 12` with `n4 + n5 + n6` at most the face cap. `max_boundary_candidate` does not drop candidates.

### d-SEAMS rings

`primitive::ringNetwork` follows Yuan and Cormack. A primitive ring of even size `2L` is a vertex plus two internally disjoint shortest paths of length `L` to its antipode. A primitive ring of odd size `2L+1` is two shortest paths of length `L` whose ends are adjacent. Each ring is emitted once, from its lowest-index vertex, and is oriented toward the lesser neighbour of that vertex. The output is sorted, so it depends on the graph and not on the order of the neighbour rows.

The candidates are then filtered by Franzblau’s test. For a cycle of length at most 7 the largest ring distance is 3, so “strictly shorter than the ring” is a walk of one or two edges (`shortcut_is_hop2`). The search reads that from the adjacency rows and does not build a bounded ball. Longer cycles still use one ball of radius `floor(k/2) − 1` per vertex, and the hop answers the pairs whose ring distance is at most 3.

Pairs that lie on a shortest path out of the source, and on an even ring the pairs that meet the antipode, cannot be the pair that fails the test (`prefix_of_shortest`, `suffix_of_shortest`, `even_antipode_dist`). They are not asked again. The ball search, which asks every pair, is kept as `primitive::ringNetworkBalls` and returns the same vectors.

### d-SEAMS cages

`cage::findBySignatures` grows connected sets of the primitive rings. A seed ring is the lowest index in its cage. A partial cage that shares two consecutive edges with a new ring, or that closes a vertex and is then offered another ring through that vertex, is cut. Up to 32 signatures share one search. The acceptance test is the same sphere condition SQQ uses: every edge twice, at least three faces in one cycle around each vertex, and `V − E + F = 2`. Ice with no pentagons makes a pentagon-bearing signature infeasible before the search branches (`ClosedSearch::feasible`).

The cage enumerator is unchanged. The ring set is unchanged, so the cages are unchanged.

## Complexity on a degree-bounded graph

Let `n` be the number of vertices, `Δ` the maximum degree (4 to 6 for a water network), and `k` the largest ring (`8` to `10` in the general bound, `6` in the ice timing). Write `L = floor(k/2)`.

| Search | Worst case | Expected, girth at least `k`, `P` shortest paths |
| --- | --- | --- |
| SQQ chordless DFS | `O(n k Δ^k)` | `O(n k Δ^k)`. No chord fires inside the ball, so the walk is the whole simple-path tree |
| SQQ GROW | exponential in the number of faces | same; one state per connected face set that fits a composition |
| d-SEAMS primitive rings, with balls | `O(n Δ^{k/2} + n B P^2 k^2)` where a vertex has at most `P` shortest paths and the ball has `B` vertices | `O(n Δ^{k/2})` when `P` is constant. The ball pass is the same order as the per-source BFS |
| d-SEAMS primitive rings, hop, `k ≤ 7` | the pairing term only; the ball pass is gone | `O(n Δ^{k/2})`, with one BFS per source instead of three |
| d-SEAMS cages | exponential in the number of faces | infeasible signatures return before branching |

`P` is the quantity that separates the two ring searches. On a worst-case graph many shortest paths share their ends and the pairing is exponential, the same order as SQQ’s path tree. On ice Ic the oxygen network is the diamond lattice: girth 6, degree 4, and each antipode of a hexagon has two shortest paths. `P` is a constant. SQQ still walks every simple path of length up to 5, which is `Θ(Δ^k)` per vertex. The primitive search walks only the breadth-first dag of depth 3, which is `Θ(Δ^{k/2})` per vertex.

So the primitive search is asymptotically faster on degree-bounded graphs of large girth, which is the water crystal. It is not asymptotically faster on a graph engineered to have exponentially many shortest paths. The two searches do not return the same rings: chordless contains primitive, and the spur hexagon is the witness.

The published timing on 110592 molecules of ice Ic is 105 ms for this search on one thread and 18.2 s for SQQ. The ratio is about 170. The asymptotic gap for `Δ = 4`, `k = 6` is `Δ^{k/2} = 8` against `Δ^k = 4096` at the branching, and the rest is the `std::set` canonicalisation, `std::map` adjacency and the recursive `std::function` in the 0.5.8 engine.

## The hop probe

`formal/Formal/Rings.lean` proves the facts the probe uses.

- `primitive_chordless`: a primitive ring has no chord.
- `hex_chordless` and `hex_not_primitive`: chordless does not imply primitive.
- `dist_le_ringDist`: an arc of the cycle is a walk, so the graph distance is at most the ring distance. Primitive means equality.
- `even_antipode_dist`: on an even primitive ring both arcs from any vertex to the opposite vertex are shortest paths. The pair enumeration is complete, and the lowest-index start emits each ring once.
- `prefix_of_shortest` and `suffix_of_shortest`: a prefix of a shortest walk is shortest, and so is the tail. Pairs on the source path and, for an even ring, pairs at the antipode are not retested.
- `dist_lt_three_iff` and `shortcut_is_hop2`: if the cycle has length at most 7, a strict shortcut is an edge or a common neighbour.
- `shortcut_within_half`: if `G` and `G'` differ by edges that touch a vertex set `C`, a primitive ring of `G` that is reachable in `G'` and has a strictly shorter pair there meets `C` within distance `k/2` in the union graph.
- `reach_of_shortcut` and `reach_zero_of_le_seven`: a shortcut of `l ≤ floor(k/2) − 1` edges meets the ring within `(l − 1)/2` steps, and that radius is 0 for `k ≤ 7`.

There is no `sorry` in `formal/Formal/Rings.lean`.

`RingUpdater` already uses the radius 0 bound for rings of size at most 7: only a vertex whose neighbour row changed can lie on a ring that appears or disappears. The naive “three-bond” ball, radius `floor(6/2) = 3`, is sound and strictly larger. The hop probe removes the matching ball refresh. For `k ≤ 7` an incremental update no longer rebuilds distance tables.

## Ice Ic, one session

Ideal ice Ic, oxygen sites of the cubic diamond cell, lattice `4 × 2.75 / sqrt(3)` Å, cutoff 3.5 Å. Every molecule has four neighbours and the second shell sits at 4.49 Å. A cell of `m` repetitions has `8 m³` molecules and `2` primitive hexagons per molecule (`12` hexagons through each molecule, each counted once). `scripts/gen_ice_ic.py` rebuilds the coordinates; the timing binary builds the same cell in memory and does not write a trajectory.

Release `-O3`, one process, warm-up discarded, six samples. `legacy` is `ringNetworkBalls` (the bounded-ball search). `hop` is `ringNetwork`. The two vectors are equal on every sample. Raw lines are in `markdown/rings-ic-timing.txt`.

| Molecules | Threads | Legacy median (ms) | Hop median (ms) | Rings |
| --- | ---: | ---: | ---: | ---: |
| 10648 | 1 | 7.746 | 5.424 | 21296 |
| 10648 | 8 | 3.628 | 3.033 | 21296 |
| 110592 | 1 | 80.847 | 58.548 | 221184 |
| 110592 | 8 | 27.160 | 23.713 | 221184 |

On this machine the ball search at 110592 molecules and one thread takes 80.8 ms. The hop search takes 58.5 ms and returns the same 221184 rings. The 105 ms figure is the earlier one-thread measurement of the ball search; the ordering is the same. Eight threads narrow the gap because the remaining work is the per-source breadth-first walk, which both searches do.

## Rings of size 8, 9 and 10

`shortcut_is_hop2` stops at length 2. The same antipodal-arc fact extends the bound: on a cycle of length at most `k`, a strict shortcut has length at most `floor(k/2) − 1`. For `k` in 8, 9 and 10 that is 3 or 4. Pairs at distance at most 3 are still read from the adjacency rows. Pairs farther out probe an open-addressed ball whose occupancy is a bitset, with a four-word Bloom filter in front of the probe. The sorted pair table is the baseline (`ringNetworkBalls`). Both return the same vectors; the random graphs in `tests/test_ring_hop.cpp` check sizes 5 through 10 against an independent breadth-first enumeration.

Ideal ice Ic has girth 6, so the primitive rings of size at most 10 are exactly the 221184 hexagons (110592 molecules) or 21296 hexagons (10648). The larger search still walks depth 4 or 5 and rejects every longer candidate. Release `-O3`, six samples after one warm-up. Raw lines are in `markdown/rings-ic-timing-v2.txt`.

| Molecules | Depth | Threads | Pair table (ms) | Bitset ball (ms) |
| --- | ---: | ---: | ---: | ---: |
| 10648 | 8 | 1 | 20.50 | 16.76 |
| 10648 | 8 | 8 | 9.07 | 9.18 |
| 110592 | 8 | 1 | 212.3 | 202.6 |
| 110592 | 8 | 8 | 81.4 | 90.6 |
| 10648 | 10 | 1 | 110.2 | 50.2 |
| 10648 | 10 | 8 | 30.7 | 17.5 |
| 110592 | 10 | 1 | 1018 | 622 |
| 110592 | 10 | 8 | 254 | 190 |

At depth 10 the bitset ball is faster on one thread and on eight. At depth 8 the one-thread run is faster and the eight-thread run is a few milliseconds slower: the ball has radius 3, and the second breadth-first fill does not always beat the packed sort once eight cores are already busy.

## SQQ 0.5.8, OVITO, cycless

Same ice, one thread, the neighbour graph or the bond table built before the timer. SQQ 0.5.8 is `find_rings` with the chordless definition and size 6. OVITO is `FindRingsModifier` after `CreateBondsModifier`. cycless is `cycles_iter` on the same undirected graph. On this crystal the chordless hexagons are the primitive ones, and every tool reported the same count.

| Molecules | d-SEAMS hop, depth 6 | OVITO | SQQ 0.5.8 | cycless |
| --- | ---: | ---: | ---: | ---: |
| 10648 | 5.47 ms | 43.8 ms | 1146 ms | 33270 ms |
| 110592 | 54.3 ms | 572 ms | 12500 ms | 353797 ms |

The 10648 foreign samples are two timed calls; the 110592 foreign samples are one timed call after a warm-up. d-SEAMS medians are six calls. OVITO at 110592 is about ten times the hop search. SQQ is about 230 times. cycless is about 6500 times.

## Cage expansion

`formal/Formal/Cage.lean` records the growth rule and what it does and does not give.

- `grows_card`: one step adds one face.
- `growth_hits_budget`: a branch that starts empty and takes at most `F` steps stays inside the face budget, so the recursion terminates.
- `walk_exit` and `connected_has_exit`: in a connected dual, every nonempty proper subset has an edge into the complement. That edge is exposed (used once by the subset, covered by the missing face). Branching on every such face therefore reaches every face of the cage.
- `canonical_parent_unique`: deleting the greatest face is the only parent, so each face set is grown once.
- `powerset_exponential`: an `F`-face dual has `2^F` subsets. A search that materialises those subsets does not have delay polynomial in `F`. The canonical parent removes duplicate cages. It does not by itself bound the dead-end branches between two solutions. No polynomial-delay enumerator follows from the exposed-edge rule alone.

## Halo of `floor(k/2)` bonds

`formal/Formal/Halo.lean` is the rank argument used by `ringHops`.

- `primitive_within_half` and `primitive_ring_subset_ball`: a primitive ring of length at most `k` lies in the closed ball of radius `floor(k/2)` about any of its vertices. A rank that owns one vertex of the ring, and keeps every vertex within that many bonds, holds the whole ring.
- `shortcut_vertex_bound` and `shortcut_within_half_radius`: a vertex of a shortcut shorter than the ring distance stays inside that same ball when the two arcs out of the owned vertex cover the cycle. The induced subgraph therefore still sees the shortcut, and a ring that is primitive on the rank is primitive in the full graph.
- `dist_along_shortcut`: the triangle inequality along a shortest walk is the step that places the shortcut vertex.

`ringHops` returns `maxDepth/2`, and each bond is at most one cutoff, so the coordinate skin of that many cutoffs is at least this bond skin.

`seams::domain::ringHaloBonds` is that bound, `kRingHaloHalf` is 2, and `ringHops` is `ringHaloBonds`. `distributedRings` passes `ringHaloBonds(maxDepth)` cutoffs to the rank partition. `tests/test_domain.cpp` checks the value for every size through 12.

There is no `sorry` in `formal/Formal/Rings.lean`, `formal/Formal/Halo.lean`, or `formal/Formal/Cage.lean`.

## The cutoff CSR

`nneigh::orthoFullCsr` writes the orthorhombic neighbour list as offsets and columns. `primitive::ringNetwork` reads that CSR. `ringsOf` uses it for an orthorhombic cloud and does not expand the rows into a vector of vectors first. A skewed cell still takes the older builder.

## Which ball, by measurement

On 110592 molecules of ice Ic the occupancy bitset lost to the sorted pair table at radius 3, which is rings of size 8 and 9: 247 ms against 205 ms on one thread, and 80 ms against 62 ms on eight, before this choice. Radius 4, rings of size 10, stayed faster as the bitset. The search now keeps the pair table for sizes 8 and 9 and the bitset for size 10. Pairs at ring distance at most 3 are still read from the adjacency rows. After that choice, size 8 on the same crystal takes 430 ms against 461 ms for the pair table alone on one thread, and 89.7 ms against 92.8 ms on eight. The ring vectors are equal.

## A shortcut through a lower index

The breadth-first walk used to ignore every vertex whose index was below the source. The comment said a shorter path through such a vertex would put the vertex on the ring. That is false: the path can be a shortcut and never a member. On the first liquid frame, 62 of the reported hexagons had graph distance 2 between an antipodal pair and ring distance 3, the length-2 path going through a lower index. The walk now uses the full graph. A ring that actually contains a vertex below the source is still dropped, because that source is not its lowest member.

With the full graph, the primitive counts are the OVITO counts. SQQ's chordless count is larger, and the extra rings are not primitive.

## Ice Ic, neighbours and rings together

Release `-O3`, one process. The d-SEAMS timer is `orthoFullCsr` plus `ringNetwork`. The OVITO timer is `CreateBondsModifier` plus `FindRingsModifier`. The SQQ timer is the cell-list graph plus `find_rings` at size 6. One warm-up, then six d-SEAMS samples and one foreign sample. Every tool returned the same ring count. Raw lines are in `markdown/rings-ic-timing-v3.txt`.

| Molecules | Threads | d-SEAMS depth 6 | d-SEAMS depth 8 | d-SEAMS depth 10 | OVITO depth 6 | SQQ depth 6 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 10648 | 1 | 31.8 ms | 45.4 ms | 161 ms | 53.8 ms | 1.30 s |
| 10648 | 8 | 3.83 ms | 7.92 ms | 23.9 ms | 51.6 ms | 1.30 s |
| 110592 | 1 | 152 ms | 452 ms | 1.88 s | 629 ms | 14.2 s |
| 110592 | 8 | 32.2 ms | 82.4 ms | 312 ms | 638 ms | 13.7 s |
| 884736 | 1 | 1.31 s | 3.61 s | 15.1 s | 5.15 s | 119 s |
| 884736 | 8 | 408 ms | 731 ms | 2.66 s | 5.05 s | 115 s |

Ring counts: 21296, 221184 and 1769472. SQQ and OVITO do not go faster from one thread to eight on this run. At 884736 molecules and eight threads OVITO takes 5.05 s and SQQ takes 115 s, against 408 ms for the depth-6 search.

## Liquid and the nucleation trajectory

The liquid is the 4096-molecule frame whose SHA-256 is `d02023892336c43cf409da8ee4b7dee2b8f30f94ee9e76a1c1d1ae99e646ea2b`, five frames. The nucleation file is figshare 11448702, 148 frames of 4096 molecules. That deposit does not contain a 13553-molecule frame, and neither does figshare project 73545. Frames 1, 50, 100 and 148 were compared, with all five liquid frames.

On every one of those frames, at depths 6 and 10, the d-SEAMS count of each size from 4 through 10 equals the OVITO histogram. Depth 8 is the same sizes, since the depth-10 histogram includes them. Sizes 3, 4 and 5 on liquid frame 1 also equal SQQ. SQQ's hexagons do not: the same frame has 5156 chordless hexagons and 2911 primitive ones. The other 2245 have no chord and fail the shortest-path test. One of them is vertices `0, 140, 2573, 426, 1756, 3818`. Vertices 0 and 426 sit three steps apart around the ring and two steps apart in the graph, and the ring has no chord. The primitive set is the one OVITO and d-SEAMS both report.

## Depth 10 under a second

On 110592 molecules the depth-10 search spent about 250 ms building a radius-4 ball and about 1.6 s pairing shortest paths. The pairing time was not the primitivity test. Most of the pairs are length 8 and 10, and on this crystal every one of them fails. A radius-2 pair table answers the questions those pairs ask: distance at most 2 is a lookup, distance 3 is a neighbour sitting at distance 2, and distance at most 4 is a vertex in both balls. Paths are kept only through vertices at or above the source, which is sound once distances are taken in the full graph. An antipode whose distance-2 keys are all within distance 2 of one another cannot close a primitive ring of size 8 or 10, so its path pairs are not built. A bipartite graph skips the odd search. The occupancy bitset of the radius-2 ball is the probe in front of that test.

`full_graph_check_necessary` in `formal/Formal/Rings.lean` is the low-index fixture. On the hexagon `1-2-3-4-5-6`, `HighOnly.dist 1 4 = 3` after vertex 0 is deleted, and `LowSpur.dist 1 4 ≤ 2` when vertex 0 is joined to 1 and 4. That is the shape of the 62 liquid hexagons. A walk that refuses to enter vertex 0 reports distance 3 and would keep the hexagon. The full graph, the brute-force primitive enumerator, and OVITO all reject it and keep the two pentagons. Chordless search keeps the hexagon as well. `shortcut_meets_radius_two` is the radius-2 meet: a strict shortcut on a cycle of length at most 10 has length at most 4, and that walk has a vertex in both radius-2 balls.

Distances are taken on the whole graph, and only then are levels below the source cleared so a ring is not extended through those vertices. The two steps are the same rule: count the shortcut, and do not put the lower vertex on the ring. The neighbour scratch is one buffer per calling thread.

Release build, one thread, eight samples after a warm-up, 8 online cores, neighbours and rings together. Median 584.7 ms, 221184 rings, the same vector as the pair-table search. The radius-2 table is 15 ms of that. Raw lines are in `markdown/rings-ic-timing-v4.txt`.

## Interface and nucleation

Serial, one thread, eight samples. The interface dump has 864 molecules. The nucleation frame has 4096. The 884736-molecule systems are a 6×6×6 tiling of the nucleation frame and a 16×8×8 tiling of the interface. Sizes 4, 5 and 6 match OVITO on every one of these. Sizes 3, 4 and 5 match SQQ on the untiled frames. SQQ's hexagons are the chordless count, which is larger. One SQQ sample on the 884736 nucleation tile did not finish within 20 minutes; the chordless search is the cost, not the neighbour build.

| System | Molecules | d-SEAMS | OVITO | SQQ |
| --- | ---: | ---: | ---: | ---: |
| Nucleation | 4096 | 27.9 ms | 33.9 ms | 561 ms |
| Nucleation tile | 884736 | 6.63 s | 9.19 s | — |
| Interface | 864 | 5.80 ms | 7.80 ms | 145 ms |
| Interface tile | 884736 | 5.91 s | 8.08 s | — |

Nucleation hexagons: d-SEAMS and OVITO 3409, SQQ 4524. Interface hexagons: d-SEAMS and OVITO 1143, SQQ 1325. The tiled counts are those numbers times the number of copies.
