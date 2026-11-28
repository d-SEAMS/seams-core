<!-- https://docs.dseams.info/explanation/ring-completion.md -->

# Ring-adjacent completion {#ring-adjacent-completion}

`ring::ringAdjacentCompletion` fills a single vacancy on a labeled
six-ring. `ring::seededCageAffiliation` has two overloads. The one
that computes the affiliation takes the completion flag as its fifth
argument, a `bool` named `ringAdjacentCompletion`, default `false`.
The overload that receives a precomputed affiliation puts the same
flag last. Cage flags are hexagonal cage (HC) and double-diamond cage
(DDC).

## The all-but-one rule {#the-all-but-one-rule}

Take the per-atom HC or DDC flags from a seeded pass and the
six-membered rings of the permissive graph. A ring whose vertices all
carry that cage label but one is a cage ring with a vacancy. The last
vertex takes the label. The walk repeats until a fixed point. The HC
flags and the DDC flags run as two walks.

A liquid six-ring that touches a nucleus has at most a few labeled
vertices. Five labeled vertices are required, so the walk cannot
enter the liquid. An empty seed stays empty.

## The flag on `seededCageAffiliation` {#the-flag-on-seededcageaffiliation}

```c++
ring::SeededAtomLabels seededCageAffiliation(
    const std::vector<std::vector<int>> &strictRings,
    const std::vector<std::vector<int>> &strictNList,
    const std::vector<std::vector<int>> &permissiveRings,
    const std::vector<std::vector<int>> &permissiveNList,
    bool ringAdjacentCompletion = false);
```

When the flag is `false` (the default), the function returns the
seeded labels. When the flag is `true`, those labels pass through
`ring::ringAdjacentCompletion` on `permissiveRings`.

`seams cages --complete` sets the flag.

On cubic ice the walk leaves both flag vectors unchanged. Clearing
one atom’s label on that crystal fills it again. An empty seed stays
empty.

## Least fixed point {#least-fixed-point}

The labeled set is the least fixed point of the all-but-one rule
above the seed. Visiting order does not change that set.

Edge-sharing is a different rule. It labels every vertex of a ring
that shares two vertices with a fully labeled ring, so two shared
vertices label the rest of that ring. Those vertices stay unlabeled
under the all-but-one closure.
