#include <domain.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>

#include <franzblau.hpp>
#include <generic.hpp>
#include <neighbours.hpp>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace seams::domain {
namespace {

// Ownership: a cube of 2^bits cells per axis in fractional coordinates, at
// least kCellsPerRank per rank so a cut falls within a small share, at most
// 2^kMaxBits per axis. Halo: cells at most halo / kHaloCells wide, so the
// halo overshoots by less than a cell, at most kMaxHaloCells per axis.
constexpr long long kCellsPerRank = 64;
constexpr int kMaxBits = 7;
constexpr double kHaloCells = 4.0;
constexpr int kMaxHaloCells = 128;
// Fewer rows or rings than this are not worth waking a team for
constexpr int kParallelMin = 4096;

//! Grid cell along one axis of a fractional coordinate
int cellAlong(double s, int n) {
  double f = s - std::floor(s);
  if (!(f >= 0.0)) {
    f = 0.0;
  }
  return std::min(n - 1, static_cast<int>(f * n));
}

//! Set every cell within @a reach cells of a set one along @a axis,
//! wrapping periodically
void dilate(std::vector<char> &mask, const int n[3], int axis, int reach) {
  const int len = n[axis];
  if (reach <= 0 || len <= 1) {
    return;
  }
  const int stride[3] = {1, n[0], n[0] * n[1]};
  const int u = (axis + 1) % 3;
  const int v = (axis + 2) % 3;
  const bool whole = 2 * reach + 1 >= len;
  std::vector<char> line(static_cast<std::size_t>(len));
  for (int a = 0; a < n[u]; a++) {
    for (int b = 0; b < n[v]; b++) {
      const int base = a * stride[u] + b * stride[v];
      bool any = false;
      for (int i = 0; i < len; i++) {
        line[i] = mask[base + i * stride[axis]];
        any = any || line[i];
      }
      if (!any) {
        continue;
      }
      for (int i = 0; i < len; i++) {
        if (whole) {
          mask[base + i * stride[axis]] = 1;
        } else if (line[i]) {
          for (int d = -reach; d <= reach; d++) {
            mask[base + ((i + d + len) % len) * stride[axis]] = 1;
          }
        }
      }
    }
  }
}

} // namespace

std::uint64_t hilbertKey(std::uint32_t x, std::uint32_t y, std::uint32_t z,
                         int bits) {
  // Skilling, AIP Conf. Proc. 707, 381 (2004): axes to transposed index
  std::uint32_t X[3] = {x, y, z};
  const std::uint32_t top = 1u << (bits - 1);
  for (std::uint32_t q = top; q > 1; q >>= 1) {
    const std::uint32_t p = q - 1;
    for (int i = 0; i < 3; i++) {
      if (X[i] & q) {
        X[0] ^= p;
      } else {
        const std::uint32_t t = (X[0] ^ X[i]) & p;
        X[0] ^= t;
        X[i] ^= t;
      }
    }
  }
  X[1] ^= X[0];
  X[2] ^= X[1];
  std::uint32_t t = 0;
  for (std::uint32_t q = top; q > 1; q >>= 1) {
    if (X[2] & q) {
      t ^= q - 1;
    }
  }
  for (int i = 0; i < 3; i++) {
    X[i] ^= t;
  }
  std::uint64_t key = 0;
  for (int b = bits - 1; b >= 0; b--) {
    for (int i = 0; i < 3; i++) {
      key = (key << 1) | ((X[i] >> b) & 1u);
    }
  }
  return key;
}

Share decompose(const Cloud &cloud, int rank, int nRanks, double halo) {
  Share share;
  const int nAtoms = static_cast<int>(cloud.pts.size());
  if (rank < 0 || rank >= std::max(nRanks, 1)) {
    return share;
  }
  const gen::FracBox b = gen::makeFracBox(cloud);
  if (nRanks <= 1 || nAtoms == 0 || !b.ok) {
    if (rank == 0) {
      share.local.resize(static_cast<std::size_t>(nAtoms));
      std::iota(share.local.begin(), share.local.end(), 0);
      share.owned.assign(static_cast<std::size_t>(nAtoms), 1);
    }
    return share;
  }

  int bits = 1;
  while (bits < kMaxBits && (1LL << (3 * bits)) < kCellsPerRank * nRanks) {
    bits++;
  }
  const int np = 1 << bits;
  const double w[3] = {b.wx, b.wy, b.wz};
  int n[3] = {1, 1, 1};
  int reach[3] = {0, 0, 0};
  for (int k = 0; k < 3; k++) {
    if (halo > 0.0) {
      n[k] = static_cast<int>(std::clamp(w[k] * kHaloCells / halo, 1.0,
                                         static_cast<double>(kMaxHaloCells)));
      reach[k] = static_cast<int>(
          std::min(std::ceil(halo * n[k] / w[k]), static_cast<double>(n[k])));
    }
  }

  // Ownership and halo cell of every atom, by fractional coordinate
  std::vector<int> partOf(static_cast<std::size_t>(nAtoms));
  std::vector<int> haloOf(static_cast<std::size_t>(nAtoms));
  std::vector<int> count(static_cast<std::size_t>(np) * np * np, 0);
  for (int i = 0; i < nAtoms; i++) {
    const auto &p = cloud.pts[static_cast<std::size_t>(i)];
    const double sz = (p.z - b.oz) / b.lz;
    const double sy = (p.y - b.oy - b.yz * sz) / b.ly;
    const double sx = (p.x - b.ox - b.xy * sy - b.xz * sz) / b.lx;
    const int c = (cellAlong(sz, np) * np + cellAlong(sy, np)) * np +
                  cellAlong(sx, np);
    partOf[static_cast<std::size_t>(i)] = c;
    haloOf[static_cast<std::size_t>(i)] =
        (cellAlong(sz, n[2]) * n[1] + cellAlong(sy, n[1])) * n[0] +
        cellAlong(sx, n[0]);
    count[static_cast<std::size_t>(c)]++;
  }

  // The keys of a cube are a permutation, so placing each cell at its key
  // orders them along the curve; cut where the running count crosses each
  // rank's equal share
  std::vector<int> along(count.size());
  for (int c = 0; c < static_cast<int>(count.size()); c++) {
    along[hilbertKey(static_cast<std::uint32_t>(c % np),
                     static_cast<std::uint32_t>((c / np) % np),
                     static_cast<std::uint32_t>(c / (np * np)), bits)] = c;
  }
  std::vector<char> owned(count.size(), 0);
  long long before = 0;
  for (const int c : along) {
    const long long here = count[static_cast<std::size_t>(c)];
    const long long owner =
        std::min<long long>(nRanks - 1, (2 * before + here) * nRanks /
                                            (2 * static_cast<long long>(nAtoms)));
    owned[static_cast<std::size_t>(c)] = owner == rank;
    before += here;
  }

  // A point within halo of an atom lies within halo/w[k] of it in fractional
  // coordinate k, which is reach[k] halo cells
  std::vector<char> near(static_cast<std::size_t>(n[0]) * n[1] * n[2], 0);
  for (int i = 0; i < nAtoms; i++) {
    if (owned[static_cast<std::size_t>(partOf[static_cast<std::size_t>(i)])]) {
      near[static_cast<std::size_t>(haloOf[static_cast<std::size_t>(i)])] = 1;
    }
  }
  for (int k = 0; k < 3; k++) {
    dilate(near, n, k, reach[k]);
  }
  for (int i = 0; i < nAtoms; i++) {
    const char mine =
        owned[static_cast<std::size_t>(partOf[static_cast<std::size_t>(i)])];
    if (mine || (halo > 0.0 &&
                 near[static_cast<std::size_t>(haloOf[static_cast<std::size_t>(i)])])) {
      share.local.push_back(i);
      share.owned.push_back(mine);
    }
  }
  return share;
}

Cloud subCloud(const Cloud &cloud, const std::vector<int> &atoms) {
  Cloud out;
  out.currentFrame = cloud.currentFrame;
  out.box = cloud.box;
  out.boxLow = cloud.boxLow;
  out.pts.reserve(atoms.size());
  for (const int i : atoms) {
    out.pts.push_back(cloud.pts[static_cast<std::size_t>(i)]);
  }
  out.nop = static_cast<int>(out.pts.size());
  out.idIndexMap.reserve(atoms.size());
  for (int k = 0; k < out.nop; k++) {
    out.idIndexMap[out.pts[static_cast<std::size_t>(k)].atomID] = k;
  }
  return out;
}

std::vector<std::vector<int>> rings(const Cloud &cloud, double cutoff,
                                    int maxDepth, int rank, int nRanks) {
  // Members of a ring lie within maxLvl hops of its source, and the search
  // reads adjacency no further out. A ball is consulted only on whether two
  // members are joined by a path shorter than their separation around the
  // ring, at most maxLvl - 1 edges, and every vertex of such a path lies
  // within half its length of one end. So no atom further than the sum of
  // the two from an owned atom bears on an owned ring. A hop spans at most
  // the cutoff.
  const int maxLvl = maxDepth / 2;
  const int hops = maxLvl + (maxLvl - 1) / 2;
  const Share share = decompose(cloud, rank, nRanks, hops * cutoff);
  if (share.local.empty()) {
    return {};
  }
  // Ascending rows: the local and the frame index orders then agree on every
  // row, and so on the order each ring is listed in
  auto nList = nneigh::getNewNeighbourListByIndex(subCloud(cloud, share.local), cutoff);
  const int nLocal = static_cast<int>(nList.size());
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nLocal >= kParallelMin && !omp_in_parallel())
#endif
  for (int k = 0; k < nLocal; k++) {
    if (nList[k].size() > 2) {
      std::sort(nList[k].begin() + 1, nList[k].end());
    }
  }
  auto out = primitive::ringNetwork(nList, maxDepth, share.owned);
  const int nRings = static_cast<int>(out.size());
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nRings >= kParallelMin && !omp_in_parallel())
#endif
  for (int r = 0; r < nRings; r++) {
    for (int &v : out[r]) {
      v = share.local[static_cast<std::size_t>(v)];
    }
  }
  return out;
}

#ifdef SEAMS_HAS_MPI
std::vector<std::vector<int>> gatherRings(const Cloud &cloud, double cutoff,
                                          int maxDepth, MPI_Comm comm) {
  int rank = 0;
  int nRanks = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nRanks);
  auto mine = rings(cloud, cutoff, maxDepth, rank, nRanks);
  if (nRanks == 1) {
    return mine;
  }
  // Each ring travels as its length followed by its members
  std::vector<int> send;
  for (const auto &ring : mine) {
    send.push_back(static_cast<int>(ring.size()));
    send.insert(send.end(), ring.begin(), ring.end());
  }
  const int sendCount = static_cast<int>(send.size());
  std::vector<int> counts(static_cast<std::size_t>(nRanks));
  std::vector<int> displs(static_cast<std::size_t>(nRanks));
  MPI_Allgather(&sendCount, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
  std::exclusive_scan(counts.begin(), counts.end(), displs.begin(), 0);
  std::vector<int> all(static_cast<std::size_t>(displs.back() + counts.back()));
  MPI_Allgatherv(send.data(), sendCount, MPI_INT, all.data(), counts.data(),
                 displs.data(), MPI_INT, comm);

  // A source's rings all come from its one owner, in the order a single
  // rank lists them, so placing them by source restores that order
  const int nAtoms = static_cast<int>(cloud.pts.size());
  std::vector<std::size_t> at;
  std::vector<std::size_t> first(static_cast<std::size_t>(nAtoms) + 1, 0);
  for (std::size_t k = 0; k < all.size(); k += static_cast<std::size_t>(all[k]) + 1) {
    at.push_back(k);
    first[static_cast<std::size_t>(all[k + 1]) + 1]++;
  }
  std::partial_sum(first.begin(), first.end(), first.begin());
  const int nRings = static_cast<int>(at.size());
  std::vector<std::size_t> slot(at.size());
  for (int r = 0; r < nRings; r++) {
    slot[r] = first[static_cast<std::size_t>(all[at[r] + 1])]++;
  }
  std::vector<std::vector<int>> out(at.size());
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nRings >= kParallelMin && !omp_in_parallel())
#endif
  for (int r = 0; r < nRings; r++) {
    const auto begin = all.begin() + static_cast<std::ptrdiff_t>(at[r]) + 1;
    out[slot[r]].assign(begin, begin + all[at[r]]);
  }
  return out;
}
#endif

} // namespace seams::domain
