/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** Order and labelling invariance of the cage predicates on one frame. The
** batch affiliation and the findHC + findDDC cage objects are evaluated on
** the ring list as enumerated and reversed, with every ring rotated by one,
** reversed, or rotated by two and reversed, and on a randomly relabelled
** copy of the frame whose rings are enumerated afresh. Each ring is keyed
** by its atom set. Reordering the ring list must give the same cages. A
** relabeling must give the same affiliation and the same cage count: the
** count is the size of a maximum non-conflicting set. The representative
** packing is chosen by the atom-index key, so the particular cages follow
** the labels.
**
** Exits 1 on any mismatch, and 2 when the trajectory cannot be read, the
** frame holds no atom of the type or the graph has no six-ring.
**
**   check_invariance TRAJ [frame] [cutoff] [atomType] [seed]
*/

#include <cage.hpp>
#include <cage_affiliation.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <numeric>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace {

using Rings = std::vector<std::vector<int>>;

std::vector<int> key(std::vector<int> ring) {
  std::sort(ring.begin(), ring.end());
  return ring;
}

std::vector<std::vector<int>> keys(const Rings &six,
                                   const std::vector<int> &ids) {
  std::vector<std::vector<int>> out;
  out.reserve(ids.size());
  for (const int i : ids) {
    out.push_back(key(six[static_cast<std::size_t>(i)]));
  }
  std::sort(out.begin(), out.end());
  return out;
}

struct Result {
  int hcRings = 0, ddcRings = 0, hcAtoms = 0, ddcAtoms = 0;
  int cages = 0, hcCages = 0, ddcCages = 0, listHC = 0, listDDC = 0;
  std::vector<int> atomHc, atomDdc;
  // Order-free forms: ring keys of each flag, list and cage
  std::vector<std::vector<int>> hcKeys, ddcKeys, listHCKeys, listDDCKeys;
  std::vector<std::pair<int, std::vector<std::vector<int>>>> cageKeys;
};

Result evaluate(const Rings &six, const Rings &nList,
                const std::vector<int> &rank) {
  Result c;
  const auto aff = ring::cageAffiliation(six, nList);
  c.atomHc.assign(nList.size(), 0);
  c.atomDdc.assign(nList.size(), 0);
  std::vector<int> hcIds, ddcIds;
  for (std::size_t i = 0; i < six.size(); i++) {
    if (aff.hc[i]) {
      hcIds.push_back(static_cast<int>(i));
    }
    if (aff.ddc[i]) {
      ddcIds.push_back(static_cast<int>(i));
    }
    for (const int a : six[i]) {
      c.atomHc[static_cast<std::size_t>(a)] |= aff.hc[i] ? 1 : 0;
      c.atomDdc[static_cast<std::size_t>(a)] |= aff.ddc[i] ? 1 : 0;
    }
  }
  c.hcRings = static_cast<int>(hcIds.size());
  c.ddcRings = static_cast<int>(ddcIds.size());
  c.hcAtoms = std::accumulate(c.atomHc.begin(), c.atomHc.end(), 0);
  c.ddcAtoms = std::accumulate(c.atomDdc.begin(), c.atomDdc.end(), 0);
  c.hcKeys = keys(six, hcIds);
  c.ddcKeys = keys(six, ddcIds);
  std::vector<ring::strucType> types(six.size(), ring::strucType::unclassified);
  std::vector<cage::Cage> cages;
  auto index = ring::buildRingSearchIndex(six, static_cast<int>(nList.size()));
  index.atomRank = rank;
  const auto listHC = ring::findHC(six, types, nList, cages, index);
  const auto listDDC = ring::findDDC(six, types, listHC, cages, index);
  c.cages = static_cast<int>(cages.size());
  for (const auto &cg : cages) {
    const bool hc = cg.type == cage::cageType::HexC;
    (hc ? c.hcCages : c.ddcCages)++;
    c.cageKeys.emplace_back(hc ? 0 : 1, keys(six, cg.rings));
  }
  std::sort(c.cageKeys.begin(), c.cageKeys.end());
  c.listHC = static_cast<int>(listHC.size());
  c.listDDC = static_cast<int>(listDDC.size());
  c.listHCKeys = keys(six, listHC);
  c.listDDCKeys = keys(six, listDDC);
  return c;
}

// Rewrites the atom-keyed fields of a result on a relabelled frame in the
// original labels
void mapBack(Result &c, const std::vector<int> &perm) {
  const std::size_t n = perm.size();
  std::vector<int> inverse(n);
  for (std::size_t a = 0; a < n; a++) {
    inverse[static_cast<std::size_t>(perm[a])] = static_cast<int>(a);
  }
  auto remap = [&](std::vector<std::vector<int>> &ks) {
    for (auto &k : ks) {
      for (int &a : k) {
        a = inverse[static_cast<std::size_t>(a)];
      }
      std::sort(k.begin(), k.end());
    }
    std::sort(ks.begin(), ks.end());
  };
  std::vector<int> hc(n), ddc(n);
  for (std::size_t a = 0; a < n; a++) {
    hc[a] = c.atomHc[static_cast<std::size_t>(perm[a])];
    ddc[a] = c.atomDdc[static_cast<std::size_t>(perm[a])];
  }
  c.atomHc = std::move(hc);
  c.atomDdc = std::move(ddc);
  remap(c.hcKeys);
  remap(c.ddcKeys);
  remap(c.listHCKeys);
  remap(c.listDDCKeys);
  for (auto &cg : c.cageKeys) {
    remap(cg.second);
  }
  std::sort(c.cageKeys.begin(), c.cageKeys.end());
}

bool sameAffiliation(const Result &a, const Result &b) {
  return a.hcRings == b.hcRings && a.ddcRings == b.ddcRings &&
         a.hcAtoms == b.hcAtoms && a.ddcAtoms == b.ddcAtoms &&
         a.atomHc == b.atomHc && a.atomDdc == b.atomDdc;
}

bool sameAll(const Result &a, const Result &b) {
  return sameAffiliation(a, b) && a.hcKeys == b.hcKeys &&
         a.ddcKeys == b.ddcKeys && a.listHCKeys == b.listHCKeys &&
         a.listDDCKeys == b.listDDCKeys && a.cageKeys == b.cageKeys;
}

void counts(const char *name, const Result &c) {
  std::printf("%-10s hcRings %d ddcRings %d hcAtoms %d ddcAtoms %d cages %d "
              "(hc %d ddc %d) listHC %d listDDC %d",
              name, c.hcRings, c.ddcRings, c.hcAtoms, c.ddcAtoms, c.cages,
              c.hcCages, c.ddcCages, c.listHC, c.listDDC);
}

} // namespace

int main(int argc, char **argv) {
  const std::string traj = argc > 1 ? argv[1] : "traj/mW_cubic.lammpstrj";
  const int frame = argc > 2 ? std::atoi(argv[2]) : 1;
  const double cutoff = argc > 3 ? std::atof(argv[3]) : 3.5;
  const int atomType = argc > 4 ? std::atoi(argv[4]) : 1;
  const unsigned seed =
      argc > 5 ? static_cast<unsigned>(std::atoi(argv[5])) : 1;

  if (!std::ifstream(traj)) {
    std::fprintf(stderr, "check_invariance: cannot read %s\n", traj.c_str());
    return 2;
  }
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrjO(traj, frame, cloud, atomType);
  if (cloud.nop == 0) {
    std::fprintf(stderr, "check_invariance: no atom of type %d in frame %d of %s\n",
                 atomType, frame, traj.c_str());
    return 2;
  }
  const auto nl = nneigh::neighListO(cutoff, cloud, atomType);
  const auto idx = nneigh::neighbourListByIndex(cloud, nl);
  Rings six;
  for (const auto &r : primitive::ringNetwork(idx, 6)) {
    if (r.size() == 6) {
      six.push_back(r);
    }
  }
  std::printf("# %s frame %d cutoff %.2f nop %d six %zu\n", traj.c_str(), frame,
              cutoff, cloud.nop, six.size());
  if (six.empty()) {
    std::fprintf(stderr, "check_invariance: no six-ring, nothing to compare\n");
    return 2;
  }
  const std::vector<int> rank = ring::positionRank(cloud);
  bool ok = true;

  Rings rot = six;
  for (auto &r : rot) {
    std::rotate(r.begin(), r.begin() + 2, r.end());
    std::reverse(r.begin(), r.end());
  }
  Rings rot1 = six;
  for (auto &r : rot1) {
    std::rotate(r.begin(), r.begin() + 1, r.end());
  }
  Rings flip = six;
  for (auto &r : flip) {
    std::reverse(r.begin(), r.end());
  }
  const std::vector<std::pair<const char *, Rings>> variants = {
      {"reversed", Rings(six.rbegin(), six.rend())},
      {"rotated", rot},
      {"rotated1", rot1},
      {"flipped", flip}};

  auto runVariants = [&](const char *order, const std::vector<int> &r,
                         const Result &ref) {
    std::printf("## claim order: %s\n", order);
    counts("enumerated", ref);
    std::printf("\n");
    for (const auto &[name, variant] : variants) {
      const Result c = evaluate(variant, idx, r);
      const bool same = sameAll(c, ref);
      counts(name, c);
      std::printf(" %s\n", same ? "same" : "DIFFERS");
      ok = ok && same;
    }
  };
  const Result ref = evaluate(six, idx, rank);
  runVariants("atom positions", rank, ref);

  // Relabel the frame and enumerate its rings afresh
  const std::size_t n = idx.size();
  std::vector<int> perm(n);
  std::iota(perm.begin(), perm.end(), 0);
  std::mt19937 rng(seed);
  std::shuffle(perm.begin(), perm.end(), rng);
  Rings pidx(n);
  for (std::size_t a = 0; a < n; a++) {
    auto &row = pidx[static_cast<std::size_t>(perm[a])];
    row.push_back(perm[a]);
    for (std::size_t k = 1; k < idx[a].size(); k++) {
      row.push_back(perm[static_cast<std::size_t>(idx[a][k])]);
    }
  }
  Rings psix;
  for (const auto &r : primitive::ringNetwork(pidx, 6)) {
    if (r.size() == 6) {
      psix.push_back(r);
    }
  }
  std::vector<int> inverse(n);
  for (std::size_t a = 0; a < n; a++) {
    inverse[static_cast<std::size_t>(perm[a])] = static_cast<int>(a);
  }
  Rings back;
  for (const auto &r : psix) {
    std::vector<int> o;
    for (const int a : r) {
      o.push_back(inverse[static_cast<std::size_t>(a)]);
    }
    back.push_back(key(o));
  }
  std::sort(back.begin(), back.end());
  Rings orig;
  for (const auto &r : six) {
    orig.push_back(key(r));
  }
  std::sort(orig.begin(), orig.end());
  const bool sameRings = back == orig;

  auto pcloud = cloud;
  for (std::size_t a = 0; a < n; a++) {
    pcloud.pts[static_cast<std::size_t>(perm[a])] = cloud.pts[a];
  }
  Result pc = evaluate(psix, pidx, ring::positionRank(pcloud));
  mapBack(pc, perm);
  const bool sameLabels = sameAffiliation(pc, ref);
  // The representative packing is the maximum set chosen by the ring key,
  // and that key reads atom indices, so a relabeling selects the packing
  // that is canonical for the new labels. The count is the invariant: it
  // is the size of a maximum non-conflicting set.
  // listHC/listDDC are the rings touched by one maximum packing. A second
  // labelling can select another packing of the same size, so those lists
  // need not match. The object counts are the invariant.
  const bool sameCounts = pc.cages == ref.cages && pc.hcCages == ref.hcCages &&
                          pc.ddcCages == ref.ddcCages;
  std::printf("relabelled six %zu, ring set %s, per-atom affiliation %s, "
              "cage count %s\n",
              psix.size(), sameRings ? "same" : "DIFFERS",
              sameLabels ? "same" : "DIFFERS", sameCounts ? "same" : "DIFFERS");
  counts("relabelled", pc);
  std::printf(" %s\n", sameCounts ? "same" : "DIFFERS");
  ok = ok && sameRings && sameLabels && sameCounts;

  runVariants("atom indices", {}, evaluate(six, idx, {}));
  std::printf("%s\n", ok ? "PASS" : "FAIL");
  return ok ? 0 : 1;
}
