<!-- https://docs.dseams.info/api/namespace_topo.md -->

# namespace topo {#namespace-topo}

## Overview {#overview}

## Detailed Documentation {#detailed-documentation}

### Typedefs {#typedefs}

Neighbour rows by index, each row leading with the atom itself, as [nneigh::neighbourListByIndex](namespace_nneigh.md#doxid-namespacenneigh-1acf46d8c4b42f15e634ff0076ac83c012) returns them.

### Global Functions {#global-functions}

Atoms within `hops` bonds of `atom`, the atom itself first, then in breadth-first order.

Weisfeiler-Lehman refinement hash of a graph given by local adjacency lists; `root` (or -1) starts in its own color; `colors` (empty or one integer per vertex, such as the atom type) seeds the initial colors.

`rounds` refinements.

Key of the rooted neighbourhood of `atom`.

`colors` is empty or one class per row (atom type, species); colored vertices never match vertices of another color.

Keys of every atom, their histogram, the ring census up to `maxRingSize`, and the frame key.

`colors` as in localKey.

Recompute local keys only for the hop-ball of `dirtyAtoms` and keep the rest of `prev`.

The ring census is taken on `rows`. A size or hop mismatch falls back to a full fingerprint.

Hex string of a 64-bit hash.

Add every distinct key of `fp` under `label`; a key already present under other labels carries all of them, sorted and joined by ‘|’.

Text form: a header line `# method M hops H colors C` (C is 0 or 1; `colours` is the deprecated alias of `colors`; a header without either word means uncolored), then `key label` lines.

Look every atom key of `fp` up in `lib`.

Throws std::invalid_argument when the methods, hop counts or colorings differ.

Match against libraries at several depths: an atom takes the label of the deepest library that holds its key, and a molecule whose wide neighbourhood is disturbed (a defect two bonds away, a jittered dense polymorph) still gets named by its inner shells.

`depth` records the hops that named each atom. Every library must share the method of this build and the coloring (`colors` empty or one per row); throws std::invalid_argument otherwise, or when two libraries share a depth.
