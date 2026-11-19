<!-- https://docs.dseams.info/explanation/topological-network.md -->

# Topological network criteria {#topological-network-criteria}

These are the criteria the classifiers use.
Worked counts are in the
[rings](../user/tutorials/rings.md),
[bulk ice](../user/tutorials/bulk-ice.md),
[order parameter](../user/tutorials/order-parameters.md) and
[confined water](../user/tutorials/confined-water.md) tutorials.
The two CHILL+ rules are compared on
[CHILL+: the published rule and the OVITO rule](chill-plus-rules.md).

## Bond correlation {#bond-correlation}

CHILL and CHILL+ look at $l = 3$.
With $N_b(i)$ neighbours,

$$
q_{3m}(i) = \frac{1}{N_b(i)} \sum_{j \in N(i)} Y_{3m}(\hat r_{ij})
$$

$$
c_{ij} = \mathrm{Re}
\frac{\sum_m q_{3m}(i)\, q_{3m}(j)^{*}}
{\sqrt{\sum_m |q_{3m}(i)|^2 \sum_m |q_{3m}(j)|^2}}
$$

A bond is staggered when $c_{ij} \le -0.8$.
The eclipsed window is the part that differs:

| rule | eclipsed | staggered | neighbours |
| --- | --- | --- | --- |
| CHILL | -0.2 \le c_{ij} \le -0.05 | c_{ij} \le -0.8 | 4 |
| CHILL+ | -0.35 \le c_{ij} \le 0.25 | c_{ij} \le -0.8 | 4 |

The code stores those edges as `chill::chillRule` `(-0.8, -0.2, -0.05, 4)` and `chill::chillPlusRule` `(-0.8, -0.35, 0.25, 4)`.
For those two rules the boundaries are inclusive.
`ovito-chill+` uses the CHILL+ numbers with open edges:
$c_{ij} < -0.8$ is staggered, and the eclipsed window is open at both ends.

The four-bond star is the phase.
Both published rules read the four nearest neighbours.
CHILL labels from that star even when a fifth oxygen lies inside the cutoff.
CHILL+ first requires the cutoff shell itself to hold exactly four neighbours, so that fifth oxygen is water.
The average in Nguyen and Molinero, equation 1 (J. Phys. Chem. B 2015, 119, 9369), is still over the four nearest,

$$
q_{lm}(i) = \frac{1}{4} \sum_{j=1}^{4} Y_{lm}(\hat r_{ij}).
$$

With four bonds the phases are:

| phase | eclipsed | staggered | further test |
| --- | --- | --- | --- |
| cubic | 0 | 4 |  |
| hexagonal | 1 | 3 |  |
| interfacial |  | 2 | a neighbour has more than two staggered bonds |
| interfacial | 0 | 3 | CHILL: a neighbour has exactly two staggered bonds. CHILL+: a neighbour has more than one |
| clathrate | 4 | 0 | CHILL+ only |
| interfacial clathrate | 3 |  | CHILL+ only |
| water |  |  | every remaining atom |

CHILL has no clathrate row.
An atom whose recorded bond count is not 4 is water under both rules.
`bop.hpp` lists the same phase names, and `chill::isInterfacial` is the neighbour test in the two interfacial rows.

`ovito-chill+` (`--rule ovito-chill+`, or the `ovito-chill+` command) is the OVITO 3.17.2 Chill+ modifier.
A crystal row still requires exactly four neighbours inside the cutoff.
$q_3$ is averaged over every neighbour in that shell, the edges are the open cuts above, and the interfacial rows do not consult a neighbour.
Two staggered bonds, or three staggered bonds and no eclipsed bond, are interfacial.

## Primitive rings {#primitive-rings}

A cycle is primitive, in the Franzblau sense, when for every pair of its vertices a shortest path in the bond graph lies along the cycle.
`primitive::ringNetwork` returns the primitive cycles whose size is at most `max_depth`.
The bond graph is an argument.
Mutual k-nearest, the union of those nominations, and a distance cutoff are three different graphs.
Quote the graph with the census.

## Cages {#cages}

A hexagonal cage is two basal six-rings and three prismatic six-rings.
A double-diamond cage is one equatorial six-ring and six peripheral rings.
`ring::cageAffiliation` evaluates those conditions for every six-ring.
`findHC` and `findDDC` claim rings as they accept them, so which rings they test depends on visiting order.
The affiliation predicates evaluate every six-ring.

The seeded variant takes a strict graph for the seeds and a permissive graph for the completions.
An atom on the permissive graph is kept only when its component contains a seed.
`--complete` then fills a six-ring that is missing one labeled vertex.
That walk is the least fixed point of the all-but-one rule, on
[Ring-adjacent completion](ring-completion.md).

## Confined ice {#confined-ice}

A prism block is two rings of the same size.
`ring::basalPrismConditions` is the strict test, and a pair that passes is perfect.
With shape matching on, a pair that fails that test is imperfect when `ring::relaxedPrismConditions` finds at least one bond between the rings and the shape match succeeds.
`ring::findPrisms` returns the rings that are blocks and the perfect and imperfect counts.

A monolayer analysis keeps primitive rings of each size and writes the coverage of the sheet.
The in-plane $g(r)$ normalises the pair histogram by the slab factor

$$
f(r) = 1 - \frac{r}{2h} \quad (r \le h), \qquad f(r) = \frac{h}{2r} \quad (r > h)
$$

where $h$ is the confined height.
That is `rdf2::normalizeRDF`.

## Order parameters {#order-parameters}

$$
q_l(i) = \sqrt{\frac{4\pi}{2l+1} \sum_{m=-l}^{l} |q_{lm}(i)|^2}
$$

$\bar q_l(i)$ replaces $q_{lm}(i)$ with the average of $q_{lm}$ over $i$ and its neighbours before the same sum.
That is the Lechner-Dellago form.
The implementation is `chill::steinhardtQl`.
It computes degrees 3, 4, 6, 8 and 12.
Other degrees return zeros.
