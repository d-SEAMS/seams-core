<!-- https://docs.dseams.info/api/namespace_nneigh.md -->

# namespace nneigh {#namespace-nneigh}

## Overview {#overview}

Functions for building neighbour lists. [More…](#details-namespacenneigh)

## Detailed Documentation {#detailed-documentation}

Functions for building neighbour lists.

This namespace contains functions that build neighbour lists (linkcell cutoff rows or k-nearest searches), saving either the atom IDs or atom indices (according to a PointCloud) in a row-ordered vector of vectors. Whether the atom IDs or atom indices (i.e. the indices of the elements in the vector pts inside the PointCloud) are saved, the neighbour lists are constructed such that the first element is the ‘central atom’, whose neighbours are being saved on that particular line. The central atom is followed by the atom IDs or indices of the nearest neighbours.

In a ‘full’ neighbour list, if 1 is a neighbour of 2 then 1 is saved in the neighbour list of 2 AND 2 is also saved in the neighbour list of 1.

In a ‘half’ neighbour list, if 1 is a neighbour of 2, then 2 is saved in the neighbour list of 1 but not vice versa.

### Changelog {#changelog}

- Amrita Goswami \[[amrita16thaug646@gmail.com](mailto:amrita16thaug646@gmail.com)\]; date modified: Nov 14, 2019
- Rohit Goswami \[[rog32@hi.is](mailto:rog32@hi.is)\]; date modified: Mar 20, 2021

#### Global Functions {#global-functions}

LAMMPS dump bound spans to restricted triclinic H (rows a, b, c).

`box[0..2]` are `xhi_bound - xlo_bound` etc. Optional `box[3..5]` are tilt `xy, xz, yz`. `boxLow` is the bound lo. Inverse of `xlo_bound = xlo + min(0, xy, xz, xy+xz)`.

Recovered restricted-triclinic lengths lx, ly, lz (H diagonal).

Longest recovered length: 0 = x, 1 = y, 2 = z.

Cartesian to fractional coordinates via dump H.

Fractional to cartesian coordinates via dump H.

Triclinic dump-cell volume, the absolute value of det(H), from dumpBoundsToH.

Bound spans Lx*Ly*Lz are not the cell volume when tilt is present.

LAMMPS edge origin and edge lengths, plus tilt, as a dump box.

`lx, ly, lz` are `xhi - xlo` and not the bound spans. A zero tilt stays three spans so the orthorhombic pair loop is unchanged. Returns false when an edge is not positive.

Crystallographic lengths and angles in degrees, stored as a dump box.

`alpha` is between b and c, `beta` between a and c, `gamma` between a and b. The basis is the reduced triclinic frame: a along x, b in the xy plane. Returns false when the cell is not a positive volume.

Map a LAMMPS dump box onto lc_cell.

Cell and optional per-frame ortho lengths for analyzeResident.

boxLow set or nBox >= 6: one dump box via lammpsBoxToLinkcell, frameLens is null. Otherwise Cell::ortho of box\[0..2\] and frameLens is box.

Cutoff rows of the cloud from linkcell’s threaded CSR search, in cloud indices.

Row `i` lists every `j != i` whose nearest image is strictly within *rcutoff*, once and ascending; *half* keeps `j > i`. An open axis of *periodic* takes no image, and a cloud with every axis open needs no box. A nonzero entry of *mask* keeps that atom; a masked atom has an empty row and is in no row. An empty mask keeps every atom. The search runs on *threads* Rayon workers, and the rows are the same at every width. *out* keeps its capacity between calls. Throws std::invalid_argument for a cutoff that is not finite and positive or a periodic axis without a box, and linkcell::Error for coordinates that are not finite.

All these functions use atom IDs and not indices.

Full neighbour list for pairs of type I and type J, from linkcell’s cutoff rows. The neighbour list does not differentiate between the types of atoms

I-J neighbour list (I==J is like-type and reuses neighListO). Unlike-type pairs come from the same linkcell rows as neighListO.

Cutoff rows by atom ID, self first, for every atom or only the atoms marked in *source*; an unmarked atom’s row is its own ID. Partners are in ascending atom ID. The row storage stays allocated on the calling thread and is refilled in place, so overlapping calls do not share it. The returned reference is valid on this thread until the next call on this thread. Same partners as neighListO for the rows that are filled.

Neighbour rows on a cell grid that covers the owned block plus its halo. Anchors are folded onto the first owned atom. Every atom, ghosts included, is shifted by a box length so it sits against that block. An axis whose empty gap is wider than the cutoff does not wrap, and the fused cutoff list runs on that window rather than on the periodic box. Returns false when the atoms fill the box or the box is not orthorhombic. Rows are atom IDs when *byId* is set, and cloud indices otherwise. Self is first, then partners in ascending order, matching neighListO. A null source marks every atom as an anchor; otherwise only marked atoms receive partners. With *everyRow* the source still places the window, and every atom in it receives partners.

Threaded cell-list neighbour rows. `rows[k]` holds the indices within `rcutoff` of `subset[k]` under the minimum image convention, ascending and without the atom itself. Atoms are binned in fractional coordinates of the recovered triclinic cell H; every axis needs at least three cells of perpendicular width `rcutoff`, so that the nearest image of every neighbour sits in one of the 27 surrounding cells. Returns false, with `rows` untouched, when the cell is too small for that. Each row is built by one thread.

When set, the cutoff lists search on one thread even if the OpenMP team is wider.

The team is left in place for the caller.

Full cutoff list as CSR of cloud indices.

Orthorhombic, or a restricted triclinic cell stored as bound spans plus xy, xz, yz.

Orthorhombic boxes (`box.size() == 3`) and restricted triclinic dump boxes (`box.size() >= 6`) both use the parallel cell sort. A null *subset* is every atom. `cols[offsets[i] .. offsets[i + 1]]` is sorted and has no self entry. A neighbour is kept when the squared minimum image is strictly below the cutoff squared. *half* keeps the pair only on the lower cloud index. The pair set does not depend on the thread count. Scratch buffers stay allocated and are refilled on the next call. Returns false when the box is not usable.

Same rows as neighListO, written into *rows*. A row keeps the capacity of the previous call, so a repeated frame does not allocate per atom.

Same rows as getNewNeighbourListByIndex, written into *rows*.

Directed k nearest neighbours as a CSR of cloud indices, no self entry.

`candidateCutoff` is the cell-size hint used by the k-nearest search. Rows are nearest first. Returns false when the cloud has no box.

Atom-ID rows, self first, from a CSR of cloud indices.

Full neighbour list for one type, from linkcell’s cutoff rows. An atom of another type keeps a row of its own ID. Empty when the cloud has no box. You can only use this for neighbour lists with one type

neighListO’s rows with the team narrowed to four threads below 4096 atoms. Empty when the cloud has no box or no atom of *typeI*. The OpenMP team is left running at that width.

Half neighbour list for one type, from linkcell’s half rows: the lower cloud index owns the pair. You can only use this for neighbour lists with one type

The following function outputs a neighbour list using indices and NOT atom IDs

Converts the neighbour list build with atom IDs into a neighbour list of atom indices, according to the pointCloud

Gets a neighbour list by index, according to a pointCloud given as the input. Assume no slices or other skullduggery

Bonded graph from the k nearest neighbours of each particle rather than a distance cutoff.

The default intersection symmetrization (mutual = true) bonds i and j only when each lists the other among its k nearest; the union alternative bonds on either nomination. In a crystal the first shell is mutual and the two coincide; on disordered packings the mutual graph is sparser, which starves accidental ring structure measured on the dense null, mutual scores zero false crystal where union reaches 2.5%. Nominations are a periodic linked-cell k-nearest search via linkcell (Allen and Tildesley): KD-trees have no minimum-image convention. candidateCutoff is only a cell-size hint: an empty hinted shell deliberately falls back to the whole periodic cell so a k-nearest query still returns the nearest available atoms. On an undistorted tetrahedral lattice with k = 4 this graph equals the first-shell cutoff graph. Rows are by atom ID with the leading self entry, like neighListO.

Same neighbours as kNearestNeighbourList, as a CSR of cloud indices with no self entry. Column order matches that list after the atom-id rows are mapped back to indices.

Packed nominations from the same search as kNearestNeighbourList. `xyz` is n by 3 and `nom` is n by k cloud indices, nearest first, with unused slots set to -1. The rows are not symmetrized.

k-nearest graph on a type set (water oxygens).

Substrate, ions and other species never appear as neighbours. An empty set is no atoms.

Mutual and union k-nearest graphs from one candidate search.

Seeded affiliation needs both; building them separately repeats the cell list.

The shell-separation certificate for the exact reduction of the k-nearest graph to the cutoff graph: returns {max_i d_k(i), min_i d_{k+1}(i)} over particles of the type.

When max_i d_k(i) <= rcutoff <= min_i d_{k+1}(i), the two graphs coincide edge for edge and every downstream graph predicate is identical. Brute force, intended for validation. A particle of the type with fewer than k others of that type returns {infinity, 0}: no cutoff gives it k neighbours.

Nearest unlike image of each typeI particle among typeJ particles.

Returns cloud-index pairs (i, j) and the MIC distance (sqrt of [gen::periodicDistSq](namespace_gen.md#doxid-namespacegen-1ac55d45ba7b50cb75763412f158ecd713)). On an ion cloud, typeI = 1 (cation vertex) and typeJ = 2 (anion vertex). This is the nearest unlike neighbour, not a coordination number.

Subset of nearestUnlike where j’s nearest typeI is i (mutual).

A contact pair is this mutual edge, not first-shell membership.

Erases memory for a vector of vectors for the neighbour list.
