//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <cage_affiliation.hpp>

#include <ring.hpp>
#include <topo_bulk.hpp>

#include <algorithm>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

//! Canonical cycle of a ring: the lexicographically least sequence over all
//! rotations and both directions, so the same ring keys identically however
//! it was enumerated
std::vector<int> canonicalKey(const std::vector<int> &ringAtoms) {
  const size_t k = ringAtoms.size();
  std::vector<int> best;
  std::vector<int> candidate(k);
  for (int dir = 0; dir < 2; dir++) {
    for (size_t start = 0; start < k; start++) {
      for (size_t m = 0; m < k; m++) {
        const size_t idx =
            (dir == 0) ? (start + m) % k : (start + k - m) % k;
        candidate[m] = ringAtoms[idx];
      }
      if (best.empty() || candidate < best) {
        best = candidate;
      }
    }
  }
  return best;
}

//! Shared per-frame machinery: the inverted index, sorted ring copies for
//! overlap tests, and cached basal/adjacency lists (one sort per ring).
struct FrameContext {
  const std::vector<std::vector<int>> &rings;
  const std::vector<std::vector<int>> &nList;
  ring::RingSearchIndex index;
  std::vector<std::vector<int>> sortedRings;
  mutable std::vector<std::vector<int>> basalCand;
  mutable std::vector<char> basalReady;
  mutable std::vector<std::vector<int>> adjRings;
  mutable std::vector<char> adjReady;
  // Scratch for findPrismatic, which appends to a list and stamps a type
  // vector neither of which affiliation reads
  std::vector<int> scratchListHC;
  std::vector<ring::strucType> scratchType;

  FrameContext(const std::vector<std::vector<int>> &rings_,
               const std::vector<std::vector<int>> &nList_)
      : rings(rings_), nList(nList_) {
    int maxAtom = 0;
    for (const auto &r : rings) {
      for (const int a : r) {
        maxAtom = std::max(maxAtom, a);
      }
    }
    index = ring::buildRingSearchIndex(rings, maxAtom + 1);
    const int nR = static_cast<int>(rings.size());
    sortedRings.resize(static_cast<std::size_t>(nR));
    for (int i = 0; i < nR; i++) {
      sortedRings[static_cast<std::size_t>(i)] =
          rings[static_cast<std::size_t>(i)];
      std::sort(sortedRings[static_cast<std::size_t>(i)].begin(),
                sortedRings[static_cast<std::size_t>(i)].end());
    }
    scratchType.assign(static_cast<std::size_t>(nR),
                       ring::strucType::unclassified);
    basalCand.resize(static_cast<std::size_t>(nR));
    adjRings.resize(static_cast<std::size_t>(nR));
    adjReady.assign(static_cast<std::size_t>(nR), 0);
    basalReady.assign(static_cast<std::size_t>(nR), 0);
  }

  [[nodiscard]] bool shareAtoms(int i, int j) const {
    const auto &a = sortedRings[i];
    const auto &b = sortedRings[j];
    size_t p = 0, q = 0;
    while (p < a.size() && q < b.size()) {
      if (a[p] == b[q]) {
        return true;
      }
      (a[p] < b[q]) ? p++ : q++;
    }
    return false;
  }

  [[nodiscard]] const std::vector<int> &ringsThroughAtom(int atom) const {
    static const std::vector<int> empty;
    if (atom < 0 ||
        atom >= static_cast<int>(index.ringsContainingAtom.size())) {
      return empty;
    }
    return index.ringsContainingAtom[atom];
  }

  [[nodiscard]] std::vector<int> collectAdjacentRings(int r) const {
    std::vector<int> out;
    for (const int atom : rings[r]) {
      for (const int s : ringsThroughAtom(atom)) {
        if (s != r) {
          out.push_back(s);
        }
      }
      if (atom >= 0 && atom < static_cast<int>(nList.size())) {
        for (size_t n = 1; n < nList[atom].size(); n++) {
          for (const int s : ringsThroughAtom(nList[atom][n])) {
            if (s != r) {
              out.push_back(s);
            }
          }
        }
      }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
  }

  //! Ordered basal-pair predicate, exactly the findHC acceptance test
  [[nodiscard]] bool basalPair(int i, int j) {
    if (i == j || shareAtoms(i, j)) {
      return false;
    }
    return ring::basalConditions(nList, rings[i], rings[j]);
  }

  [[nodiscard]] const std::vector<int> &adjacentRings(int r) const {
    static const std::vector<int> empty;
    if (r < 0 || r >= static_cast<int>(adjRings.size())) {
      return empty;
    }
    if (!adjReady[static_cast<std::size_t>(r)]) {
      adjRings[static_cast<std::size_t>(r)] = collectAdjacentRings(r);
      adjReady[static_cast<std::size_t>(r)] = 1;
    }
    return adjRings[static_cast<std::size_t>(r)];
  }

  //! Candidate second-basal rings for i in the findHC gathering order:
  //! rings holding a bond-neighbour of i's first or second atom
  [[nodiscard]] std::vector<int> collectBasalCandidates(int i) const {
    std::vector<int> out;
    for (int slot = 0; slot < 2 &&
                       slot < static_cast<int>(rings[i].size());
         slot++) {
      const int anchor = rings[i][slot];
      if (anchor < 0 || anchor >= static_cast<int>(nList.size())) {
        continue;
      }
      for (size_t n = 1; n < nList[anchor].size(); n++) {
        for (const int s : ringsThroughAtom(nList[anchor][n])) {
          if (s != i) {
            out.push_back(s);
          }
        }
      }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
  }

  [[nodiscard]] const std::vector<int> &basalCandidates(int i) const {
    static const std::vector<int> empty;
    if (i < 0 || i >= static_cast<int>(basalCand.size())) {
      return empty;
    }
    if (!basalReady[static_cast<std::size_t>(i)]) {
      basalCand[static_cast<std::size_t>(i)] = collectBasalCandidates(i);
      basalReady[static_cast<std::size_t>(i)] = 1;
    }
    return basalCand[static_cast<std::size_t>(i)];
  }

  //! Prismatic rings of the ordered passing pair (i, j)
  std::vector<int> prismaticOf(int i, int j) {
    std::vector<int> prismatic;
    scratchListHC.clear();
    (void)ring::findPrismatic(rings, scratchListHC, scratchType, i, j,
                              prismatic, index);
    return prismatic;
  }

  //! The equatorial-conditions test of findDDC without the claim skip; on a
  //! pass, fills the six unique peripheral ring indices
  bool equatorialPass(int i, bool hcAffiliated,
                      std::vector<int> &peripherals) {
    if (hcAffiliated) {
      return false;
    }
    peripherals.clear();
    if (!ring::conditionOneDDC(rings, peripherals, i, index)) {
      return false;
    }
    if (!ring::conditionTwoDDC(rings, peripherals, i, index)) {
      return false;
    }
    if (!ring::conditionThreeDDC(rings, peripherals)) {
      return false;
    }
    std::sort(peripherals.begin(), peripherals.end());
    peripherals.erase(std::unique(peripherals.begin(), peripherals.end()),
                      peripherals.end());
    return peripherals.size() == 6;
  }
};

//! HC affiliation of one ring from its own perspective: basal of a passing
//! ordered pair in either direction, or prismatic of a passing pair drawn
//! from its adjacency neighbourhood (any pair with a prismatic ring k has
//! both members inside A(k))
bool hcAffiliatedOf(FrameContext &ctx, int r,
                    const std::vector<int> &adjacency) {
  for (const int s : adjacency) {
    if (ctx.basalPair(r, s) || ctx.basalPair(s, r)) {
      return true;
    }
  }
  for (const int i : adjacency) {
    if (i == r) {
      continue;
    }
    const std::vector<int> &candidates = ctx.basalCandidates(i);
    for (const int j : candidates) {
      if (j == r) {
        continue;
      }
      // Both members of a pair with r prismatic lie in A(r)
      if (!std::binary_search(adjacency.begin(), adjacency.end(), j)) {
        continue;
      }
      if (!ctx.basalPair(i, j)) {
        continue;
      }
      const std::vector<int> prismatic = ctx.prismaticOf(i, j);
      if (std::find(prismatic.begin(), prismatic.end(), r) !=
          prismatic.end()) {
        return true;
      }
    }
  }
  return false;
}

} // namespace

ring::CageAffiliation
ring::cageAffiliation(const std::vector<std::vector<int>> &rings,
                      const std::vector<std::vector<int>> &nList) {
  CageAffiliation result;
  const size_t nRings = rings.size();
  result.hc.assign(nRings, false);
  result.ddc.assign(nRings, false);
  if (nRings == 0) {
    return result;
  }
  FrameContext ctx(rings, nList);

  // HC affiliation by a global sweep over ordered pairs: every passing pair
  // marks its basal rings and its prismatic rings. Sweeping every i covers
  // both directions of every pair once.
  for (size_t i = 0; i < nRings; i++) {
    const std::vector<int> &candidates = ctx.basalCandidates(static_cast<int>(i));
    for (const int j : candidates) {
      if (!ctx.basalPair(static_cast<int>(i), j)) {
        continue;
      }
      result.hc[i] = true;
      result.hc[j] = true;
      for (const int k : ctx.prismaticOf(static_cast<int>(i), j)) {
        result.hc[k] = true;
      }
    }
  }

  // DDC affiliation: every ring passing the equatorial conditions marks
  // itself and its six peripherals
  std::vector<int> peripherals;
  for (size_t i = 0; i < nRings; i++) {
    if (ctx.equatorialPass(static_cast<int>(i), result.hc[i], peripherals)) {
      result.ddc[i] = true;
      for (const int p : peripherals) {
        result.ddc[p] = true;
      }
    }
  }
  return result;
}

// ---------------------------------------------------------------------------

struct ring::AffiliationUpdater::Impl {
  bool primed = false;
  std::vector<std::vector<int>> prevRings;
  std::vector<std::vector<int>> prevSortedRows; // nList rows, sorted
  std::map<std::vector<int>, std::pair<bool, bool>>
      storedByKey; // canonical key -> (hc, ddc)
  std::vector<std::vector<int>> keys; // canonical key per previous ring
  CageAffiliation current;
  int lastReclassified = 0;

  void primeFull(const std::vector<std::vector<int>> &rings,
                 const std::vector<std::vector<int>> &nList) {
    current = ring::cageAffiliation(rings, nList);
    lastReclassified = static_cast<int>(rings.size());
    remember(rings, nList);
    primed = true;
  }

  void remember(const std::vector<std::vector<int>> &rings,
                const std::vector<std::vector<int>> &nList) {
    prevRings = rings;
    prevSortedRows.resize(nList.size());
    for (size_t a = 0; a < nList.size(); a++) {
      prevSortedRows[a] = nList[a];
      std::sort(prevSortedRows[a].begin(), prevSortedRows[a].end());
    }
    keys.resize(rings.size());
    storedByKey.clear();
    for (size_t i = 0; i < rings.size(); i++) {
      keys[i] = canonicalKey(rings[i]);
      storedByKey[keys[i]] = {current.hc[i], current.ddc[i]};
    }
  }
};

ring::AffiliationUpdater::AffiliationUpdater() : impl_(new Impl) {}
ring::AffiliationUpdater::~AffiliationUpdater() = default;
ring::AffiliationUpdater::AffiliationUpdater(AffiliationUpdater &&) noexcept =
    default;
ring::AffiliationUpdater &
ring::AffiliationUpdater::operator=(AffiliationUpdater &&) noexcept = default;

int ring::AffiliationUpdater::lastReclassified() const {
  return impl_->lastReclassified;
}

const ring::CageAffiliation &
ring::AffiliationUpdater::update(const std::vector<std::vector<int>> &rings,
                                 const std::vector<std::vector<int>> &nList) {
  Impl &st = *impl_;
  if (!st.primed || st.prevSortedRows.size() != nList.size()) {
    st.primeFull(rings, nList);
    return st.current;
  }

  // ------------------------------------------------------------------ diff
  // Seed atoms: atoms of added or removed rings, plus atoms whose neighbour
  // row changed. Unchanged rows are identical in both frames, so expanding
  // the seeds through the new rows reaches everything the old graph reached.
  std::vector<std::vector<int>> newKeys(rings.size());
  std::unordered_set<int> seedAtoms;
  {
    std::map<std::vector<int>, int> newByKey;
    for (size_t i = 0; i < rings.size(); i++) {
      newKeys[i] = canonicalKey(rings[i]);
      newByKey.emplace(newKeys[i], static_cast<int>(i));
    }
    for (size_t i = 0; i < rings.size(); i++) { // added rings
      if (st.storedByKey.find(newKeys[i]) == st.storedByKey.end()) {
        seedAtoms.insert(rings[i].begin(), rings[i].end());
      }
    }
    for (size_t i = 0; i < st.prevRings.size(); i++) { // removed rings
      if (newByKey.find(st.keys[i]) == newByKey.end()) {
        seedAtoms.insert(st.prevRings[i].begin(), st.prevRings[i].end());
      }
    }
    std::vector<int> row;
    for (size_t a = 0; a < nList.size(); a++) { // changed neighbour rows
      row = nList[a];
      std::sort(row.begin(), row.end());
      if (row != st.prevSortedRows[a]) {
        seedAtoms.insert(static_cast<int>(a));
      }
    }
  }

  if (seedAtoms.empty()) {
    // Same topology; only the ring order may differ. Carry stored answers.
    st.current.hc.assign(rings.size(), false);
    st.current.ddc.assign(rings.size(), false);
    for (size_t i = 0; i < rings.size(); i++) {
      const auto it = st.storedByKey.find(newKeys[i]);
      st.current.hc[i] = it->second.first;
      st.current.ddc[i] = it->second.second;
    }
    st.lastReclassified = 0;
    st.remember(rings, nList);
    return st.current;
  }

  FrameContext ctx(rings, nList);

  // Dirty closure: two rounds of (atoms -> one bond-hop -> rings -> atoms),
  // matching the A(A(r)) dependence of the DDC predicate
  auto expandOnce = [&](const std::unordered_set<int> &atoms) {
    std::unordered_set<int> dirty;
    for (const int a : atoms) {
      for (const int r : ctx.ringsThroughAtom(a)) {
        dirty.insert(r);
      }
      if (a >= 0 && a < static_cast<int>(nList.size())) {
        for (size_t n = 1; n < nList[a].size(); n++) {
          for (const int r : ctx.ringsThroughAtom(nList[a][n])) {
            dirty.insert(r);
          }
        }
      }
    }
    return dirty;
  };
  std::unordered_set<int> dirty1 = expandOnce(seedAtoms);
  std::unordered_set<int> atoms2 = seedAtoms;
  for (const int r : dirty1) {
    atoms2.insert(rings[r].begin(), rings[r].end());
  }
  std::unordered_set<int> dirty = expandOnce(atoms2);
  dirty.insert(dirty1.begin(), dirty1.end());

  // ----------------------------------------------------------- reclassify
  st.current.hc.assign(rings.size(), false);
  st.current.ddc.assign(rings.size(), false);
  for (size_t i = 0; i < rings.size(); i++) {
    if (dirty.count(static_cast<int>(i))) {
      continue;
    }
    const auto it = st.storedByKey.find(newKeys[i]);
    // A clean ring existed in the previous frame: it is not added, and its
    // locality closure is untouched
    st.current.hc[i] = it->second.first;
    st.current.ddc[i] = it->second.second;
  }

  // Per-ring recomputation on the new frame for the dirty closure. Cache the
  // helper answers shared between neighbouring dirty rings.
  std::unordered_map<int, bool> hcCache;
  auto hcOf = [&](int r) {
    const auto it = hcCache.find(r);
    if (it != hcCache.end()) {
      return it->second;
    }
    const bool value = hcAffiliatedOf(ctx, r, ctx.adjacentRings(r));
    hcCache.emplace(r, value);
    return value;
  };
  std::unordered_map<int, std::pair<bool, std::vector<int>>> eqCache;
  auto eqOf = [&](int r) -> const std::pair<bool, std::vector<int>> & {
    auto it = eqCache.find(r);
    if (it == eqCache.end()) {
      std::vector<int> peripherals;
      const bool pass = ctx.equatorialPass(r, hcOf(r), peripherals);
      it = eqCache.emplace(r, std::make_pair(pass, std::move(peripherals)))
               .first;
    }
    return it->second;
  };

  for (const int r : dirty) {
    st.current.hc[r] = hcOf(r);
    // DDC: r passes itself, or a ring sharing an atom with r passes and
    // lists r among its peripherals
    bool ddc = eqOf(r).first;
    if (!ddc) {
      for (const int atom : rings[r]) {
        for (const int i : ctx.ringsThroughAtom(atom)) {
          if (i == r) {
            continue;
          }
          const auto &eq = eqOf(i);
          if (eq.first && std::find(eq.second.begin(), eq.second.end(), r) !=
                              eq.second.end()) {
            ddc = true;
            break;
          }
        }
        if (ddc) {
          break;
        }
      }
    }
    st.current.ddc[r] = ddc;
  }

  st.lastReclassified = static_cast<int>(dirty.size());
  st.remember(rings, nList);
  return st.current;
}

ring::SeededAtomLabels ring::seededCageAffiliation(
    const std::vector<std::vector<int>> &strictRings,
    const std::vector<std::vector<int>> &strictNList,
    const std::vector<std::vector<int>> &permissiveRings,
    const std::vector<std::vector<int>> &permissiveNList) {
  const int nAtoms = static_cast<int>(permissiveNList.size());
  SeededAtomLabels out;
  out.hc.assign(nAtoms, false);
  out.ddc.assign(nAtoms, false);
  if (nAtoms == 0) {
    return out;
  }

  const auto toAtoms = [nAtoms](const std::vector<std::vector<int>> &rings,
                                const CageAffiliation &affiliation) {
    std::pair<std::vector<bool>, std::vector<bool>> flags{
        std::vector<bool>(nAtoms, false), std::vector<bool>(nAtoms, false)};
    for (size_t i = 0; i < rings.size(); i++) {
      for (const int atom : rings[i]) {
        if (atom < 0 || atom >= nAtoms) {
          continue;
        }
        if (affiliation.hc[i]) {
          flags.first[atom] = true;
        }
        if (affiliation.ddc[i]) {
          flags.second[atom] = true;
        }
      }
    }
    return flags;
  };

  const auto [hcS, ddcS] =
      toAtoms(strictRings, cageAffiliation(strictRings, strictNList));
  const auto [hcP, ddcP] = toAtoms(
      permissiveRings, cageAffiliation(permissiveRings, permissiveNList));

  // One flood per cage type: an HC seed must not keep a DDC-only atom
  // that happens to sit in the same permissive H-bond component.
  auto flood = [&](const std::vector<bool> &seeds,
                   const std::vector<bool> &affiliated) {
    std::vector<bool> kept(nAtoms, false);
    std::vector<int> stack;
    for (int a = 0; a < nAtoms; a++) {
      if (!seeds[a] || kept[a]) {
        continue;
      }
      stack.push_back(a);
      kept[a] = true;
      while (!stack.empty()) {
        const int u = stack.back();
        stack.pop_back();
        if (u < 0 || static_cast<size_t>(u) >= permissiveNList.size()) {
          continue;
        }
        for (size_t m = 1; m < permissiveNList[u].size(); m++) {
          const int v = permissiveNList[u][m];
          if (v >= 0 && v < nAtoms && affiliated[v] && !kept[v]) {
            kept[v] = true;
            stack.push_back(v);
          }
        }
      }
    }
    return kept;
  };

  std::vector<bool> affiliatedHC(nAtoms, false);
  std::vector<bool> affiliatedDDC(nAtoms, false);
  for (int a = 0; a < nAtoms; a++) {
    affiliatedHC[a] = hcS[a] || hcP[a];
    affiliatedDDC[a] = ddcS[a] || ddcP[a];
  }
  const auto keptHC = flood(hcS, affiliatedHC);
  const auto keptDDC = flood(ddcS, affiliatedDDC);

  for (int a = 0; a < nAtoms; a++) {
    if (keptHC[a]) {
      out.hc[a] = hcS[a] || hcP[a];
    }
    if (keptDDC[a]) {
      out.ddc[a] = ddcS[a] || ddcP[a];
    }
  }
  return out;
}
