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
 *  The result is identical to cageAffiliation on the new frame. When
 *  more than half the rings sit in the dirty closure, the update is
 *  the batch classifier.
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
 * @param[in] strictRings Six-membered rings of the strict graph.
 * @param[in] strictNList Strict graph, by index with leading self entries.
 * @param[in] permissiveRings Six-membered rings of the permissive graph.
 * @param[in] permissiveNList Permissive graph, same conventions.
 */
[[nodiscard]] SeededAtomLabels seededCageAffiliation(
    const std::vector<std::vector<int>> &strictRings,
    const std::vector<std::vector<int>> &strictNList,
    const std::vector<std::vector<int>> &permissiveRings,
    const std::vector<std::vector<int>> &permissiveNList);

} // namespace ring

#endif // SEAMS_CAGE_AFFILIATION_H_
