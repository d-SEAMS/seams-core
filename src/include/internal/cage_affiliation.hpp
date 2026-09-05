//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#ifndef SEAMS_CAGE_AFFILIATION_H_
#define SEAMS_CAGE_AFFILIATION_H_

#include <memory>
#include <vector>

/** @file cage_affiliation.hpp
 *  @brief Order-free per-ring cage classification with an exact incremental
 *  update.
 *
 *  The greedy cage assembly in findHC/findDDC claims rings as it accepts
 *  cages, so which rings are *tested* depends on visiting order. The
 *  affiliation predicates here are claim-free restatements of the same
 *  geometric conditions, evaluated for every ring independently:
 *
 *  - A ring is HC-affiliated when it is a basal ring of some ordered pair
 *    passing the basal conditions, or a prismatic ring of such a pair.
 *  - A ring is DDC-affiliated when it passes the equatorial conditions
 *    (excluding HC-affiliated rings from candidacy, as the published scheme
 *    does), or is one of the six peripheral rings of a ring that passes.
 *
 *  Locality theorem. Let two six-membered rings be *adjacent* when one
 *  contains an atom within one bond-hop of an atom of the other (a symmetric
 *  relation, written A(r) for the set adjacent to r). Every quantity the
 *  predicates read is confined to a bounded adjacency neighbourhood:
 *
 *  1. An ordered basal pair P(i,j) requires an atom of j to neighbour the
 *     first or second atom of i, so j is in A(i); evaluating P reads only
 *     rings i, j and the neighbour rows of their atoms.
 *  2. A prismatic ring k of a pair (i,j) shares a triplet with i and three
 *     further atoms with j, so both i and j lie in A(k).
 *  3. The equatorial conditions for ring i read only rings sharing an atom
 *     with i -- a subset of A(i) -- plus the HC affiliation of i itself.
 *  4. A peripheral ring r of an equator i shares an atom with i, so i lies
 *     in A(r).
 *
 *  Hence HC affiliation of r is a function of the rings in A(r) and the
 *  neighbour rows of their atoms, and DDC affiliation of r is a function of
 *  the rings in A(A(r)) and their rows. Between two frames, affiliation can
 *  therefore change only for rings whose second adjacency neighbourhood
 *  touches an added ring, a removed ring, or an atom whose neighbour row
 *  changed. The incremental updater recomputes exactly that closure and
 *  carries every other ring's stored answer, so its output equals a full
 *  recomputation by construction.
 */

namespace ring {

/** @struct CageAffiliation
 * @brief Per-ring affiliation flags, indexed like the input ring vector.
 */
struct CageAffiliation {
  std::vector<bool> hc;   //! Basal or prismatic of a passing basal pair
  std::vector<bool> ddc;  //! Equatorial pass, or peripheral of one
};

/** Six-rings that are ice-I stacking planes.
 *
 *  A basal ring of a passing HC pair is a hexagonal stacking plane.
 *  A DDC equatorial ring (HC-affiliated rings excluded, as in findDDC)
 *  is a cubic stacking plane. Prismatic and peripheral rings connect
 *  planes and do not vote. The two plane sets are disjoint because an
 *  HC-affiliated ring is not an equatorial candidate.
 */
struct StackingPlanes {
  std::vector<bool> basal;       //! HC basal of a passing pair
  std::vector<bool> equatorial;  //! DDC equatorial pass
};

[[nodiscard]] StackingPlanes stackingPlanes(
    const std::vector<std::vector<int>> &rings,
    const std::vector<std::vector<int>> &nList);

//! Batch computation of the affiliation predicates for six-membered rings.
//! nList is by index with the leading self entry, as neighbourListByIndex
//! produces.
[[nodiscard]] CageAffiliation cageAffiliation(
    const std::vector<std::vector<int>> &rings,
    const std::vector<std::vector<int>> &nList);

/** @class AffiliationUpdater
 * @brief Exact incremental affiliation across frames.
 * @details Diffs the six-ring set (by canonical ring key) and the neighbour
 *  rows against the previous frame, recomputes affiliation inside the
 *  locality closure of the changes, and carries stored answers elsewhere.
 *  The result is identical to cageAffiliation on the new frame. The
 *  neighbour graph should come from SkinNeighborList (default: mutual
 *  four-nearest, the TUM v2 graph).
 */
class AffiliationUpdater {
public:
  AffiliationUpdater();
  ~AffiliationUpdater();
  AffiliationUpdater(AffiliationUpdater &&) noexcept;
  AffiliationUpdater &operator=(AffiliationUpdater &&) noexcept;
  AffiliationUpdater(const AffiliationUpdater &) = delete;
  AffiliationUpdater &operator=(const AffiliationUpdater &) = delete;

  //! Classify the frame, reusing the previous frame's answers outside the
  //! locality closure of what changed. The reference stays valid until the
  //! next update call.
  const CageAffiliation &update(const std::vector<std::vector<int>> &rings,
                                const std::vector<std::vector<int>> &nList);

  //! Rings reclassified by the last update (the dirty closure); the first
  //! call reports every ring.
  [[nodiscard]] int lastReclassified() const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

/** @struct SeededAtomLabels
 * @brief Per-atom cage flags from seeded (hysteresis) affiliation.
 */
struct SeededAtomLabels {
  std::vector<bool> hc;  //! Atom belongs to an accepted HC-affiliated ring
  std::vector<bool> ddc; //! Atom belongs to an accepted DDC-affiliated ring
};

/** Seeded affiliation over two graphs on the same atoms: the strict graph
 *  (typically the mutual k-nearest bonds) supplies seeds, the permissive
 *  supergraph (typically the union bonds) supplies completions, and a
 *  permissively affiliated atom is accepted only when its bonded component
 *  of affiliated atoms contains a seed. Specificity on structureless input
 *  is structural rather than statistical: when the strict pass affiliates
 *  nothing, nothing is accepted regardless of what the permissive graph
 *  builds. Where both graphs label an atom, the strict labels win.
 *
 *  The fifth argument is the ring-adjacent completion flag. The default
 *  is false. When the flag is true the accepted labels pass through
 *  ringAdjacentCompletion() on the permissive six-rings. HC and DDC
 *  complete on separate flag vectors. seams cages --complete sets the
 *  flag. tests/walk_compare leaves the flag false.
 * @param[in] strictRings Six-membered rings of the strict graph.
 * @param[in] strictNList Strict graph, by index with leading self entries.
 * @param[in] permissiveRings Six-membered rings of the permissive graph.
 * @param[in] permissiveNList Permissive graph, same conventions.
 * @param[in] ringAdjacentCompletion When true, call
 *  ringAdjacentCompletion() on the accepted labels. Default false.
 */
[[nodiscard]] SeededAtomLabels seededCageAffiliation(
    const std::vector<std::vector<int>> &strictRings,
    const std::vector<std::vector<int>> &strictNList,
    const std::vector<std::vector<int>> &permissiveRings,
    const std::vector<std::vector<int>> &permissiveNList,
    bool ringAdjacentCompletion = false);

/** Ring completion of seeded labels.
 *
 *  A permissive six-ring whose vertices all carry a cage label but one
 *  is a cage ring with a vacancy. The last vertex takes that label.
 *  The walk repeats until a fixed point. HC and DDC complete on
 *  separate flag vectors.
 *
 *  A liquid ring that only touches a nucleus has at most a few labelled
 *  vertices, so the all-but-one rule cannot walk into the liquid. An
 *  empty seed stays empty: no ring then has five labelled vertices.
 *
 *  The result is the least fixed point above the seed. Visiting order
 *  does not change the labelled set. lean/DseamsProofs/Completion.lean
 *  states that claim. Catch2 covers the C++ identity on a crystal, the
 *  refill of a single vacancy, and the structural zero.
 * @param[in] labels Seeded per-atom flags.
 * @param[in] permissiveRings Six-membered rings of the permissive graph.
 */
[[nodiscard]] SeededAtomLabels
ringAdjacentCompletion(const SeededAtomLabels &labels,
                       const std::vector<std::vector<int>> &permissiveRings);

} // namespace ring

#endif // SEAMS_CAGE_AFFILIATION_H_
