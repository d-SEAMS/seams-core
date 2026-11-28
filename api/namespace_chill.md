<!-- https://docs.dseams.info/api/namespace_chill.md -->

# namespace chill {#namespace-chill}

## Overview {#overview}

CHILL and CHILL+ structure classification. [More…](#details-namespacechill)

## Detailed Documentation {#detailed-documentation}

CHILL and CHILL+ structure classification.

This namespace contains functions that are used in the CHILL/CHILL+ classification scheme, as well as a yodaCloud struct to hold all the information.

Both the [CHILL](https://pubs.rsc.org/en/content/articlehtml/2010/cp/b919724a) and [CHILL+](https://pubs.acs.org/doi/abs/10.1021/jp510289t) methods are based on the local bond order parameter method (ten Wolde et al., J. Chem. Phys. 104, 9932, 1996, DOI 10.1063/1.471721), developed for the identification of crystal nuclei in Lennard-Jones systems. The local order parameter method is based on [order parameters](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.47.1297) introduced by Steinhardt et al.

The local environment of the water molecules is classified using an algorithm based on spherical harmonics, which is independent of the specific crystal structure and does not require the definition of a reference frame. The local order around each water molecule $i$ is described by a local orientational bond order parameter $q_{lm}(i)$, with $2l+1$ complex components:

$$
q_{lm}(i) = @frac{1}{N_b(i)} \Sigma_{j=1}^{N_b(i)} Y_{lm}(r_{ij})
$$

Here, $N_b(i)=4$ is the number of nearest neighbours for water molecule $i$, $l$ is a free integer parameter, $m$ is an integer such that $m \in [-l,l]$. The functions $Y_{lm}(r_{ij})$ are the spherical harmonics and $r_{ij}$ is the distance vector between molecule $i$ and each one of its four nearest neighbours $j$.

The alignment of the orientation of the local structures is measured by the normalized dot product of the local orientational bond order parameter, $q_l(i)$, between each pair of neighbor molecules $i$ and $j$, given by:

$$
a(i,j) = \frac{\sum_{m=-l}^{l} q_{lm}(i)\, q_{lm}^*(j)}
{|q_l(i)|\,|q_l(j)|}
$$

where, $q_{lm}^*$ is the complex conjugate of $q_{lm}$.

Depending on the values of $a(i,j)$, the bond between each pair of molecules $i$ and $j$ is classified as being either staggered or eclipsed.

Since each molecule $i$ in deeply supercooled water has four nearest neighbours, the type of the resultant four bonds is used to identify phases according to the CHILL and CHILL+ algorithms.

The CHILL algorithm can classify water molecules as belonging to the cubic, hexagonal, interfacial or liquid (amorphous) phase:

Here, E and S refer to eclipsed and staggered bonds. The neighbours column contains the number of nearest neighbours.

The CHILL+ algorithm, which is a modified version of CHILL, additionally identifies clathrate and interfacial clathrate phases. The criteria are enumerated in the table below:

Although both the CHILL and CHILL+ classification schemes take into account the local order of the environment of each particle, both schemes output a per-particle classification.

### Changelog {#changelog}

- Amrita Goswami \[[amrita16thaug646@gmail.com](mailto:amrita16thaug646@gmail.com)\]; date modified: Sept 19, 2019

#### Global Functions {#global-functions}

Staggered, eclipsed, or out of range. *closed* selects inclusive edges (c <= staggeredMax, eclipsedMin <= c <= eclipsedMax), the Nguyen and Molinero rule and the 2.13 behavior. Open edges match OVITO 3.17.2 (c < staggeredMax, eclipsedMin < c < eclipsedMax). The flag is an argument so [BondClassifier](struct_chill_BondClassifier.md#doxid-structchill-1-1-bond-classifier) stays three doubles and one int. A caller that passes the struct by value does not carry a fifth field.

The CHILL rule set (Moore et al., PCCP 12, 4124, 2010): strict staggered bound, narrow eclipsed band, four nearest neighbours

The CHILL+ rule set (Nguyen and Molinero, JPCB 119, 9369, 2015): wider eclipsed band for interfacial and clathrate recognition

Look up a registered rule set by name. “CHILL” and “CHILL+” are always present; registerBondClassifier adds more. Throws std::out_of_range for unknown names.

Register (or replace) a named rule set at runtime, making the material’s thresholds available to every front end by name

Names of every registered rule set.

Compute and classify the bond correlations c_ij under an arbitrary rule set, filling yCloud.pts\[i\].c_ij. getCorrel and getCorrelPlus are this engine under their canonical water rules.

Same classification as the four-argument form. *closed* false uses the open OVITO edges. The four-argument form stays the 2.13 symbol and keeps inclusive edges.

CHILL/CHILL+ bonds from a CSR of cloud indices.

Bonds and CHILL+ labels from a cutoff CSR. Returns the cubic count. *keepBonds* false still writes iceType and skips the per-atom bond vectors. The labels are the same either way. The vectors are what a caller that reads c_ij needs.

OVITO Chill+ from a cutoff CSR. q_lm is the sum over every neighbour in the row, not the four closest. Edges are open and the Table 1 neighbour test is off, which is the Chill+ modifier. A crystal row still needs exactly four neighbours. Bond vectors are not stored. Returns the cubic count.

Published CHILL labels from a cutoff CSR. Correlations are the four nearest neighbours in the list. Bond vectors are not stored. Returns the cubic count. Thresholds are [chillRule()](#doxid-namespacechill-1a1bc0f302f0567e60061fd3f82a472b81).

c_ij for atoms *which* selects, from a bond-order vector that already holds every neighbour the bond reads. Atoms outside *which* are left untouched, including a ghost’s bond list.

Function for getting the bond order correlations $c_{ij}$ (or $a_{ij}$ in some treatments) according to the CHILL algorithm.

**Parameters:**

| yCloud | The output molSys::PointCloud |
| --- | --- |
| nList | The row-ordered neighbour list, by ID. The first element of each row is the particle ID, followed by the IDs of the neighbours |
| isSlice | This decides whether there is a slice or not |
| coordinationNumber | How many nearest neighbours the bond sum runs over. The default of four is the coordination CHILL is defined and validated against (tetrahedral water); other values reuse the c_ij machinery for differently coordinated systems, and a non-positive value gives each atom its own count from its neighbour-list row. The ice-classification tables in getIceType remain four-bond rules |

Function that classifies every particle’s [molSys::atom_state_type](enum_molSys_atom_state_type.md#doxid-namespacemol-sys-1a13ffe3a39a0816627553f195e102abb5) ice type, according to the CHILL algorithm.

Does not print out the information.

**Parameters:**

| yCloud | The output molSys::PointCloud |
| --- | --- |
| isSlice | This decides whether there is a slice or not |
| nList | Row-ordered neighbour list by atom ID |

Function that classifies every particle’s [molSys::atom_state_type](enum_molSys_atom_state_type.md#doxid-namespacemol-sys-1a13ffe3a39a0816627553f195e102abb5) ice type, according to the CHILL algorithm.

**Parameters:**

| yCloud | The output molSys::PointCloud |
| --- | --- |
| nList | Row-ordered neighbour list by atom ID |
| path | Path to the output directory to which ice types are written out to |
| firstFrame | First frame to be analyzed |
| isSlice | This decides whether there is a slice or not |
| outputFileName | Name of the output file, to which the ice types will be written out. |

Gets c_ij and then classifies bond types according to the CHILL+ algorithm. coordinationNumber as in getCorrel: four is the validated CHILL+ water scheme, non-positive keeps each atom’s whole neighbour row. This is the 2.13 symbol: inclusive edges (closedBounds true).

Same classification. closedBounds false uses the open OVITO edges.

Published CHILL+ on the cutoff shell. The pair rows are read inside this call; the shared neighbour-list builder is not modified. Returns false when that search is not available.

One cutoff list and one bond pass for local q4, local q6, their Lechner-Dellago averages, and the published CHILL+ labels. The four q arrays match steinhardtQlCutoff with averaging on. Returns false when the atoms are not one type or the cutoff list is not available, and then writes nothing.

Cubic pre-filter from discovery/pysr_prefilter.py. The selected expression is cmean^4. It writes cubic and returns true only when every atom has a cutoff shell of four and a negative mean correlation whose fourth power sits above the highest non-cubic score in that search. That forces every bond below -0.84, which the published rule labels staggered, so the atom is cubic. Otherwise it returns false and leaves the cloud unchanged. The default published call does not use it.

Classifies each atom according to the CHILL+ algorithm. This is the 2.13 symbol: the published neighbour test stays on.

neighbourTest false is the OVITO interfacial rule.

CHILL+ ice types on the cloud. Does not write a file. This is the 2.13 symbol: the published neighbour test stays on.

neighbourTest false is the OVITO interfacial rule (two staggered bonds, or three staggered and no eclipsed bond, with no neighbour lookup). Either rule still requires a cutoff shell of four.

Set iceType to unclassified for atoms whose LAMMPS type is not in types. Empty types leaves the cloud unchanged. Used after a water-type k-NN so substrate and ions are not scored as water.

q6 can distinguish between water and ice. Use this for the largest ice cluster

‘Test’ condition for classifying hexagonal ice using averaged q6 and q3 Checks water According to [https://pubs.rsc.org/en/content/articlehtml/2011/cp/c1cp22167a](https://pubs.rsc.org/en/content/articlehtml/2011/cp/c1cp22167a) Gets c_ij and then classifies bond types according to the CHILL+ algorithm

Prints out the iceType for a particular frame onto the terminal.

Interfacial ice. Moore CHILL uses a neighbour with exactly two staggered bonds; Nguyen CHILL+ uses a neighbour with more than one. Pass chillPlus=true for the latter.

Finds the number of staggered bonds for a given atom of index jatom.

Local and neighbour-averaged Steinhardt parameters of degree 3, 4, 6, 8, or 12. Local ql is Steinhardt, Nelson and Ronchetti; qlBar is the Lechner-Dellago average of q_lm over the particle and its neighbours. l=12 is host-only (associated Legendre, or sphericart when linked) for the Zeron hydrate pair with l=3; the device Ylm path has no l=12. Other degrees return zeros. The compute path flattens the neighbour list to CSR so the same kernel can run on the host (OpenMP), across MPI ranks (atom split plus Allgatherv of q_lm), or under OpenMP target offload when the build provides a device.

Local $q_l$ on the cutoff list.

Neighbours are orthoFullCsr, not a four-nearest rebuild, and the kernel reads that CSR directly. Lechner-Dellago $\bar{q}_l$ is filled only when *average* is true. Otherwise `qlBar` stays zero. *invariants* fills `wl` with the normalized Steinhardt $\hat{w}_l$ for l = 4 or 6.

Local q4 and q6 from one cutoff list. Same values as two steinhardtQlCutoff calls with averaging off. Returns false when the shared shell cannot be built.

q4 and q6 from one cutoff list. The bond vectors are wrapped once. Each degree then uses the same harmonic as steinhardtQlCutoff, so the values match two separate calls. Returns false when that list is not available; the outputs are left untouched.

Mutual k-nearest Steinhardt parameters. The coordinate buffer built for the nomination is the one the harmonics read, and the neighbour rows stay in one packed array. l is 3, 4 or 6; other degrees return the result of building the mutual graph and calling steinhardtQl.

Whether steinhardtQl splits atoms over the ranks of an initialised MPI_COMM_WORLD (true unless set). A driver that gives each rank whole frames of its own turns it off.

Width of one atom’s interleaved q_lm: real, imaginary, m = -l..l.

q_lm for atoms *which* selects. Other atoms are left zero.

ql and qlBar for atoms *which* selects, reading q_lm of their neighbours.

Host CHILL+ (classifyBonds + getIceTypePlusNoPrint).

Flattened CHILL+.

Target device when compiled with offload.

Neighbours in one cutoff row, ignoring the leading atom id and any id that is not in the cloud. This is the Table 1 shell count.

`chill+`, `ovito-chill+`, a YAML document with `rule:` or `chill-plus:`, or a Lua table fragment `rule = "..."`. An empty string is the published rule. Anything else throws.

Read `spec` as a rule name. If `spec` names a readable file, the file is scanned for a `rule` or `chill-plus` key. Comments are skipped. The file is not executed as Lua or YAML.

Published or OVITO CHILL+ on one cloud. The neighbour list is the cutoff shell the rule counts.

Overlay the k nearest neighbours of each particle onto FCC, HCP, BCC and SC shells. Uses IRA when linked; otherwise Horn on a distance-sorted correspondence (high-symmetry shells only).

Same overlay, reading a CSR of cloud indices.

Bartok SOAP power spectrum for one particle (nMax radial Gaussians, spherical harmonics through lMax). Length is nMax*nMax*(lMax+1).

Bartok SOAP power spectrum for every particle. Length is nop; each row is nMax*nMax*(lMax+1).

\[q4, q6, q8\] from the Voronoi-weighted Steinhardt path (Mickel).

\[q4, q6, q8\] for every particle. One tessellation per order (l = 4, 6, 8).

Facet neighbours and area weights for every particle, with an a posteriori exactness certificate.

candidateCutoff seeds the bisector search; the cutoff is enlarged automatically until the certificate holds or the growth cap is reached.

Certificate. The cell of particle i clipped against every candidate within cutoff c is the true Voronoi cell whenever every vertex of the clipped cell lies within c/2 of i: a particle k beyond the cutoff has abs(d_k) > c, its bisector plane sits at distance abs(d_k)/2 > c/2 from i, and a half-space whose boundary is farther than c/2 cannot cut a region contained in the ball of radius c/2. An open or under-clipped cell keeps vertices on the seeding square at distance ~c and fails the certificate, which is exactly the failure the enlargement retries.

Steinhardt parameters of degree orderL (3, 4, 6 or 8) with the bond sum weighted by Voronoi facet areas. qlBar averages the weighted q_lm over the particle and its facet neighbours, after Lechner and Dellago.
