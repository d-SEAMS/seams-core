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
