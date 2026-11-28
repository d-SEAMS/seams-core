<!-- https://docs.dseams.info/reference/glossary.md -->

# Glossary {#glossary}

| term | meaning |
| --- | --- |
| libyodaLib | C++ engine library in this repository |
| seams | Argum CLI in this repository |
| yoda | Compiled Python module in PydSEAMSlib (_core / cyoda aliases) |
| dseams | Lua table from require("dseams") in yodaStruct |
| CHILL+ | Four-neighbour staggered/eclipsed star classifier |
| CHILL | Original four-neighbour rule (seams chill) |
| HC | Hexagonal cage (ice Ih, 12 waters) |
| DDC | Double-diamond cage (ice Ic, 14 waters) |
| ice score | Per-atom HC / DDC membership (seams cages) |
| seeded | Mutual four-nearest seeds, union completions |
| ringAdjacentCompletion | All-but-one fill of a labeled six-ring; fifth argument of seededCageAffiliation, default false; CLI --complete |
| walk_compare | tests/walk_compare: per-frame CHILL+ and cage counts, with the largest cluster of each score |
| nop | Number of particles in the analysed cloud |
| TUM | Topological unit matching |
| PointCloud | Positions, types, box, and per-atom labels |
| linkcell | The neighbour backend, always built. lc_pairs_within_csr builds every cutoff list. nominateByLinkcell is the k-NN for knn and seeded graphs, and for the device TUM batch |
