#include <domain.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <numeric>
#include <unordered_map>

#include <franzblau.hpp>
#include <generic.hpp>
#include <neighbours.hpp>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace seams::domain {
namespace {

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

// Factor nRanks into (nx, ny, nz). The product of owned/covered along each
// axis is the owned fraction of a uniform box; the most compact block wins
// a tie, which is the cut with the least surface.
struct BlockGrid {
  int n[3];
};

BlockGrid chooseGrid(int nRanks, const double w[3], double halo) {
  BlockGrid best{{nRanks, 1, 1}};
  double bestScore = -1.0;
  double bestMin = -1.0;
  for (int nx = 1; nx <= nRanks; nx++) {
    if (nRanks % nx != 0) {
      continue;
    }
    const int nyz = nRanks / nx;
    for (int ny = 1; ny <= nyz; ny++) {
      if (nyz % ny != 0) {
        continue;
      }
      const int nz = nyz / ny;
      const int n[3] = {nx, ny, nz};
      double score = 1.0;
      double minEdge = 1.0e300;
      for (int k = 0; k < 3; k++) {
        const double edge = w[k] > 0.0 ? w[k] : 1.0;
        const double thick = edge / n[k];
        minEdge = std::min(minEdge, thick);
        const double covered = std::min(edge, thick + 2.0 * std::max(halo, 0.0));
        score *= thick / covered;
      }
      if (score > bestScore + 1.0e-12 ||
          (std::abs(score - bestScore) <= 1.0e-12 && minEdge > bestMin)) {
        bestScore = score;
        bestMin = minEdge;
        best = {{nx, ny, nz}};
      }
    }
  }
  return best;
}

void blockAt(int rank, const int n[3], int out[3]) {
  out[0] = rank % n[0];
  out[1] = (rank / n[0]) % n[1];
  out[2] = rank / (n[0] * n[1]);
}

// Periodic gap from fractional coordinate @a f to the half-open interval
// [a, b) on the unit circle.
double fracGap(double f, double a, double b) {
  const double width = b - a;
  if (width >= 1.0 - 1.0e-15) {
    return 0.0;
  }
  double x = f - a;
  x -= std::floor(x);
  if (x <= width) {
    return 0.0;
  }
  return std::min(x - width, 1.0 - x);
}

double blockDist2(const double s[3], const int blk[3], const int n[3],
                  const double w[3]) {
  double sum = 0.0;
  for (int k = 0; k < 3; k++) {
    const double a = static_cast<double>(blk[k]) / n[k];
    const double b = static_cast<double>(blk[k] + 1) / n[k];
    const double gap = fracGap(s[k], a, b) * w[k];
    sum += gap * gap;
  }
  return sum;
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
      share.owner.assign(static_cast<std::size_t>(nAtoms), 0);
    }
    return share;
  }

  const double w[3] = {b.wx, b.wy, b.wz};
  const BlockGrid grid = chooseGrid(nRanks, w, halo);
  int mineBlk[3];
  blockAt(rank, grid.n, mineBlk);
  const double halo2 = halo > 0.0 ? halo * halo : -1.0;

  for (int i = 0; i < nAtoms; i++) {
    const auto &p = cloud.pts[static_cast<std::size_t>(i)];
    const double sz = (p.z - b.oz) / b.lz;
    const double sy = (p.y - b.oy - b.yz * sz) / b.ly;
    const double sx = (p.x - b.ox - b.xy * sy - b.xz * sz) / b.lx;
    const double s[3] = {sx, sy, sz};
    const int ix = cellAlong(sx, grid.n[0]);
    const int iy = cellAlong(sy, grid.n[1]);
    const int iz = cellAlong(sz, grid.n[2]);
    const int owner = ix + grid.n[0] * (iy + grid.n[1] * iz);
    const bool owns = owner == rank;
    if (owns || (halo2 >= 0.0 && blockDist2(s, mineBlk, grid.n, w) <= halo2)) {
      share.local.push_back(i);
      share.owned.push_back(owns ? 1 : 0);
      share.owner.push_back(owner);
    }
  }
  return share;
}

void markSources(const Cloud &cloud, const Share &share, double inner,
                 std::vector<char> &source) {
  source.assign(share.local.size(), 0);
  if (share.local.empty()) {
    return;
  }
  if (inner < 0.0) {
    for (std::size_t k = 0; k < share.local.size(); k++) {
      source[k] = share.owned[k];
    }
    return;
  }
  int rank = 0;
  int nRanks = 1;
  for (std::size_t k = 0; k < share.local.size(); k++) {
    nRanks = std::max(nRanks, share.owner[k] + 1);
    if (share.owned[k]) {
      rank = share.owner[k];
    }
  }
  const gen::FracBox box = gen::makeFracBox(cloud);
  const double w[3] = {box.wx, box.wy, box.wz};
  const BlockGrid grid = chooseGrid(nRanks, w, inner);
  int blk[3];
  blockAt(rank, grid.n, blk);
  const double inner2 = inner * inner;
  for (std::size_t k = 0; k < share.local.size(); k++) {
    if (share.owned[k]) {
      source[k] = 1;
      continue;
    }
    const auto &p = cloud.pts[k];
    const double sz = (p.z - box.oz) / box.lz;
    const double sy = (p.y - box.oy - box.yz * sz) / box.ly;
    const double sx = (p.x - box.ox - box.xy * sy - box.xz * sz) / box.lx;
    const double s[3] = {sx, sy, sz};
    if (blockDist2(s, blk, grid.n, w) <= inner2) {
      source[k] = 1;
    }
  }
}

std::vector<int> haloRanks(const Share &share) {
  std::vector<int> ranks;
  ranks.reserve(share.local.size());
  for (std::size_t k = 0; k < share.local.size(); k++) {
    if (!share.owned.empty() && !share.owned[k]) {
      ranks.push_back(share.owner[k]);
    }
  }
  std::sort(ranks.begin(), ranks.end());
  ranks.erase(std::unique(ranks.begin(), ranks.end()), ranks.end());
  return ranks;
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

namespace {
//! Rings on a cloud already limited to one share, by local index.
std::vector<std::vector<int>> ringsOf(const Cloud &local, double cutoff,
                                      int maxDepth, const std::vector<char> &owned) {
  if (local.nop == 0) {
    return {};
  }
  // Orthorhombic frames take the cutoff CSR directly. The ring search
  // reads those offsets and columns; it does not expand them into a
  // vector of rows first.
  {
    std::vector<int> offsets;
    std::vector<int> cols;
    if (nneigh::orthoFullCsr(local, nullptr, cutoff, offsets, cols)) {
      return primitive::ringNetwork(offsets, cols, maxDepth, owned);
    }
  }
  static std::vector<std::vector<int>> nList;
  static std::vector<std::vector<int>> idRows;
  static std::vector<int> idToIndex;
  if (nneigh::blockWindowRows(local, cutoff, nullptr, idRows)) {
    if (static_cast<int>(nList.size()) != local.nop) {
      nList.assign(static_cast<std::size_t>(local.nop), {});
    } else {
      for (auto &row : nList) {
        row.clear();
      }
    }
    int maxId = 0;
    for (int k = 0; k < local.nop; k++) {
      maxId = std::max(maxId, local.pts[static_cast<std::size_t>(k)].atomID);
    }
    idToIndex.assign(static_cast<std::size_t>(maxId) + 1, -1);
    for (int k = 0; k < local.nop; k++) {
      const int id = local.pts[static_cast<std::size_t>(k)].atomID;
      if (id >= 0) {
        idToIndex[static_cast<std::size_t>(id)] = k;
      }
    }
    for (int k = 0; k < local.nop; k++) {
      auto &row = nList[static_cast<std::size_t>(k)];
      row.push_back(k);
      for (std::size_t j = 1; j < idRows[static_cast<std::size_t>(k)].size(); j++) {
        const int id = idRows[static_cast<std::size_t>(k)][j];
        const int idx = (id >= 0 && id <= maxId) ? idToIndex[static_cast<std::size_t>(id)] : -1;
        if (idx >= 0) {
          row.push_back(idx);
        }
      }
    }
  } else {
    nList = nneigh::getNewNeighbourListByIndex(local, cutoff);
  }
  const int nLocal = static_cast<int>(nList.size());
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nLocal >= kParallelMin && !omp_in_parallel())
#endif
  for (int k = 0; k < nLocal; k++) {
    if (nList[k].size() > 2) {
      std::sort(nList[k].begin() + 1, nList[k].end());
    }
  }
  return primitive::ringNetwork(nList, maxDepth, owned);
}

} // namespace

Share shareWithin(const Cloud &cloud, int rank, int nRanks, double gridHalo,
                  double reach) {
  Share share;
  const int nAtoms = static_cast<int>(cloud.pts.size());
  if (rank < 0 || rank >= std::max(nRanks, 1) || nAtoms == 0) {
    return share;
  }
  const gen::FracBox b = gen::makeFracBox(cloud);
  if (!b.ok) {
    return share;
  }
  const double w[3] = {b.wx, b.wy, b.wz};
  const BlockGrid grid = chooseGrid(nRanks, w, gridHalo);
  int mineBlk[3];
  blockAt(rank, grid.n, mineBlk);
  const double reach2 = reach > 0.0 ? reach * reach : -1.0;
  for (int i = 0; i < nAtoms; i++) {
    const auto &p = cloud.pts[static_cast<std::size_t>(i)];
    const double sz = (p.z - b.oz) / b.lz;
    const double sy = (p.y - b.oy - b.yz * sz) / b.ly;
    const double sx = (p.x - b.ox - b.xy * sy - b.xz * sz) / b.lx;
    const double s[3] = {sx, sy, sz};
    const int ix = cellAlong(sx, grid.n[0]);
    const int iy = cellAlong(sy, grid.n[1]);
    const int iz = cellAlong(sz, grid.n[2]);
    const int owner = ix + grid.n[0] * (iy + grid.n[1] * iz);
    const bool owns = owner == rank;
    if (owns || (reach2 >= 0.0 && blockDist2(s, mineBlk, grid.n, w) <= reach2)) {
      share.local.push_back(i);
      share.owned.push_back(owns ? 1 : 0);
      share.owner.push_back(owner);
    }
  }
  return share;
}

namespace {
//! Largest ring size through each local atom, 1 through none
void localTypes(const std::vector<std::vector<int>> &ringList,
                std::vector<int> &type) {
  for (const auto &ring : ringList) {
    const int size = static_cast<int>(ring.size());
    for (const int v : ring) {
      type[static_cast<std::size_t>(v)] =
          std::max(type[static_cast<std::size_t>(v)], size);
    }
  }
}
} // namespace

std::vector<std::vector<int>> ownedRings(const Cloud &local, double cutoff,
                                         int maxDepth,
                                         const std::vector<char> &owned) {
  return ringsOf(local, cutoff, maxDepth, owned);
}

std::vector<std::vector<int>>
ringsFromRows(const std::vector<int> &frame,
              const std::vector<std::vector<int>> &neigh,
              const std::vector<int> &sources, int maxDepth) {
  std::unordered_map<int, int> rowAt;
  rowAt.reserve(frame.size() * 2 + 1);
  const int nRows = std::min(static_cast<int>(frame.size()),
                             static_cast<int>(neigh.size()));
  for (int i = 0; i < nRows; i++) {
    rowAt[frame[static_cast<std::size_t>(i)]] = i;
  }
  // Vertices are the frames whose rows were supplied. A neighbour past
  // that set is an edge the coordinate halo would not have contained.
  std::vector<int> verts = frame;
  verts.insert(verts.end(), sources.begin(), sources.end());
  std::sort(verts.begin(), verts.end());
  verts.erase(std::unique(verts.begin(), verts.end()), verts.end());
  std::unordered_map<int, int> localOf;
  localOf.reserve(verts.size() * 2 + 1);
  for (int i = 0; i < static_cast<int>(verts.size()); i++) {
    localOf[verts[static_cast<std::size_t>(i)]] = i;
  }
  const int n = static_cast<int>(verts.size());
  std::vector<std::vector<int>> nList(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    auto &dst = nList[static_cast<std::size_t>(i)];
    dst.push_back(i);
    const auto it = rowAt.find(verts[static_cast<std::size_t>(i)]);
    if (it == rowAt.end()) {
      continue;
    }
    for (const int f : neigh[static_cast<std::size_t>(it->second)]) {
      const auto jt = localOf.find(f);
      if (jt == localOf.end() || jt->second == i) {
        continue;
      }
      dst.push_back(jt->second);
    }
    if (dst.size() > 2) {
      std::sort(dst.begin() + 1, dst.end());
    }
  }
  std::vector<char> mask(static_cast<std::size_t>(n), 0);
  for (const int s : sources) {
    const auto it = localOf.find(s);
    if (it != localOf.end()) {
      mask[static_cast<std::size_t>(it->second)] = 1;
    }
  }
  auto rings = primitive::ringNetwork(nList, maxDepth, mask);
  for (auto &ring : rings) {
    for (int &v : ring) {
      v = verts[static_cast<std::size_t>(v)];
    }
  }
  return rings;
}

std::vector<std::vector<int>> rings(const Cloud &cloud, double cutoff,
                                    int maxDepth, int rank, int nRanks) {
  if (nRanks <= 1 || cloud.nop == 0) {
    std::vector<char> all(static_cast<std::size_t>(std::max(cloud.nop, 0)), 1);
    return ringsOf(cloud, cutoff, maxDepth, all);
  }
  // Rows out to ringRowHops are what the owners exchange. Building them
  // from the whole frame here is the same list an owner computes inside
  // its own one-cutoff halo.
  const Share near =
      shareWithin(cloud, rank, nRanks, cutoff, ringRowHops(maxDepth) * cutoff);
  if (near.local.empty()) {
    return {};
  }
  const auto full = nneigh::getNewNeighbourListByIndex(cloud, cutoff);
  std::vector<int> frame;
  std::vector<std::vector<int>> neigh;
  std::vector<int> sources;
  frame.reserve(near.local.size());
  neigh.reserve(near.local.size());
  for (std::size_t k = 0; k < near.local.size(); k++) {
    const int f = near.local[k];
    frame.push_back(f);
    std::vector<int> row;
    if (f >= 0 && static_cast<std::size_t>(f) < full.size() && !full[f].empty()) {
      row.assign(full[static_cast<std::size_t>(f)].begin() + 1,
                 full[static_cast<std::size_t>(f)].end());
    }
    neigh.push_back(std::move(row));
    if (near.owned[k]) {
      sources.push_back(f);
    }
  }
  return ringsFromRows(frame, neigh, sources, maxDepth);
}

RingSummary summarize(const std::vector<std::vector<int>> &ringList,
                      int nAtoms, int maxDepth) {
  RingSummary s;
  s.sizeCount.assign(static_cast<std::size_t>(std::max(maxDepth, 0)) + 1, 0);
  s.atomType.assign(static_cast<std::size_t>(nAtoms), 1);
  for (const auto &ring : ringList) {
    if (ring.size() < s.sizeCount.size()) {
      s.sizeCount[ring.size()]++;
    }
  }
  localTypes(ringList, s.atomType);
  return s;
}

#ifdef SEAMS_HAS_MPI
namespace {
struct GhostAtom {
  int index = 0;
  int atomID = 0;
  int type = 0;
  int molID = 0;
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

std::vector<int> symmetricNeighbors(std::vector<int> nbrs, int rank, int nRanks,
                                   MPI_Comm comm) {
  const int mine = static_cast<int>(nbrs.size());
  std::vector<int> degs(static_cast<std::size_t>(nRanks), 0);
  MPI_Allgather(&mine, 1, MPI_INT, degs.data(), 1, MPI_INT, comm);
  std::vector<int> off(static_cast<std::size_t>(nRanks), 0);
  std::exclusive_scan(degs.begin(), degs.end(), off.begin(), 0);
  const int total = nRanks == 0 ? 0 : off.back() + degs.back();
  std::vector<int> listed(static_cast<std::size_t>(std::max(total, 0)));
  MPI_Allgatherv(nbrs.empty() ? nullptr : nbrs.data(), mine, MPI_INT, listed.data(),
                 degs.data(), off.data(), MPI_INT, comm);
  std::vector<char> seen(static_cast<std::size_t>(nRanks), 0);
  for (const int n : nbrs) {
    if (n >= 0 && n < nRanks) {
      seen[static_cast<std::size_t>(n)] = 1;
    }
  }
  for (int other = 0; other < nRanks; other++) {
    if (other == rank) {
      continue;
    }
    const int base = off[static_cast<std::size_t>(other)];
    const int nEdge = degs[static_cast<std::size_t>(other)];
    for (int e = 0; e < nEdge; e++) {
      if (listed[static_cast<std::size_t>(base + e)] == rank &&
          !seen[static_cast<std::size_t>(other)]) {
        nbrs.push_back(other);
        seen[static_cast<std::size_t>(other)] = 1;
      }
    }
  }
  std::sort(nbrs.begin(), nbrs.end());
  nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
  return nbrs;
}
} // namespace

struct BoundaryPlan {
  int nRanks = 0;
  int degree = 0;
  MPI_Comm graph = MPI_COMM_NULL;
  std::vector<int> sendCount;
  std::vector<int> sendLocal;
  std::vector<int> recvCount;
  std::vector<int> recvLocal;
};

BoundaryPlan gPlan;
double gValueCommMs = 0.0;
bool gInflight = false;
double *gValues = nullptr;
int gStride = 0;
std::vector<double> gSendPack;
std::vector<double> gRecvPack;
std::vector<int> gSendDoubles;
std::vector<int> gSendDoubleDispl;
std::vector<int> gRecvDoubles;
std::vector<int> gRecvDoubleDispl;
MPI_Request gValueReq = MPI_REQUEST_NULL;

Work &gatheredWork() {
  static Work work;
  return work;
}

Work &gatherWork(const Cloud &cloud, double halo, MPI_Comm comm) {
  Work &work = gatheredWork();
  int rank = 0;
  int nRanks = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nRanks);
  if (nRanks <= 1) {
    gPlan.degree = 0;
    gPlan.nRanks = nRanks;
    if (work.cloud.nop != cloud.nop || work.share.local.size() != cloud.pts.size()) {
      work.share = decompose(cloud, rank, nRanks, halo);
      work.cloud = subCloud(cloud, work.share.local);
    } else {
      for (int i = 0; i < cloud.nop; i++) {
        work.cloud.pts[static_cast<std::size_t>(i)].x = cloud.pts[static_cast<std::size_t>(i)].x;
        work.cloud.pts[static_cast<std::size_t>(i)].y = cloud.pts[static_cast<std::size_t>(i)].y;
        work.cloud.pts[static_cast<std::size_t>(i)].z = cloud.pts[static_cast<std::size_t>(i)].z;
      }
    }
    return work;
  }

  // The neighbour relation depends only on the blocks and the halo, so the
  // graph is built once per process and reused for later frames.
  struct GraphCache {
    int nRanks = -1;
    int rank = -1;
    int nop = -1;
    double halo = -1.0;
    double haloFraction = 0.0;
    Share share;
    std::vector<int> nbrs;
    std::vector<std::vector<int>> sendAt;
    MPI_Comm graph = MPI_COMM_NULL;
  };
  static GraphCache cache;
  const bool reuse = cache.nop == cloud.nop && cache.nRanks == nRanks &&
                     cache.rank == rank && cache.halo == halo;
  if (!reuse) {
    work.share = decompose(cloud, rank, nRanks, halo);
    int nOwned = 0;
    for (const char flag : work.share.owned) {
      nOwned += flag ? 1 : 0;
    }
    const int nLocal = static_cast<int>(work.share.local.size());
    work.haloFraction =
        nLocal > 0 ? static_cast<double>(nLocal - nOwned) / nLocal : 0.0;
    if (cache.graph != MPI_COMM_NULL) {
      MPI_Comm_free(&cache.graph);
    }
    cache.nRanks = nRanks;
    cache.rank = rank;
    cache.nop = cloud.nop;
    cache.halo = halo;
    cache.share = work.share;
    cache.haloFraction = work.haloFraction;
    cache.nbrs = symmetricNeighbors(haloRanks(work.share), rank, nRanks, comm);
    const int cachedDegree = static_cast<int>(cache.nbrs.size());
    MPI_Dist_graph_create_adjacent(
        comm, cachedDegree, cache.nbrs.empty() ? nullptr : cache.nbrs.data(),
        MPI_UNWEIGHTED, cachedDegree,
        cache.nbrs.empty() ? nullptr : cache.nbrs.data(), MPI_UNWEIGHTED,
        MPI_INFO_NULL, 0, &cache.graph);
    const gen::FracBox box = gen::makeFracBox(cloud);
    const double w[3] = {box.wx, box.wy, box.wz};
    const BlockGrid grid = chooseGrid(nRanks, w, halo);
    const double halo2 = halo > 0.0 ? halo * halo : -1.0;
    cache.sendAt.assign(static_cast<std::size_t>(cachedDegree), {});
    for (std::size_t k = 0; k < work.share.local.size(); k++) {
      if (!work.share.owned[k]) {
        continue;
      }
      const auto &p = cloud.pts[static_cast<std::size_t>(work.share.local[k])];
      const double sz = (p.z - box.oz) / box.lz;
      const double sy = (p.y - box.oy - box.yz * sz) / box.ly;
      const double sx = (p.x - box.ox - box.xy * sy - box.xz * sz) / box.lx;
      const double s[3] = {sx, sy, sz};
      for (int sidx = 0; sidx < cachedDegree; sidx++) {
        int blk[3];
        blockAt(cache.nbrs[static_cast<std::size_t>(sidx)], grid.n, blk);
        if (halo2 >= 0.0 && blockDist2(s, blk, grid.n, w) <= halo2) {
          cache.sendAt[static_cast<std::size_t>(sidx)].push_back(static_cast<int>(k));
        }
      }
    }
    gPlan.nRanks = nRanks;
    gPlan.degree = cachedDegree;
    gPlan.graph = cache.graph;
    gPlan.sendCount.resize(static_cast<std::size_t>(cachedDegree));
    gPlan.sendLocal.clear();
    gPlan.recvLocal.clear();
    gPlan.recvCount.clear();
    for (int sidx = 0; sidx < cachedDegree; sidx++) {
      gPlan.sendCount[static_cast<std::size_t>(sidx)] =
          static_cast<int>(cache.sendAt[static_cast<std::size_t>(sidx)].size());
      const auto &at = cache.sendAt[static_cast<std::size_t>(sidx)];
      gPlan.sendLocal.insert(gPlan.sendLocal.end(), at.begin(), at.end());
    }
  } else {
    work.share = cache.share;
    work.haloFraction = cache.haloFraction;
  }
  const std::vector<int> &nbrs = cache.nbrs;
  const int degree = static_cast<int>(nbrs.size());
  MPI_Comm graph = cache.graph;
  std::vector<std::vector<GhostAtom>> bags(static_cast<std::size_t>(degree));
  for (int sidx = 0; sidx < degree; sidx++) {
    auto &bag = bags[static_cast<std::size_t>(sidx)];
    bag.reserve(cache.sendAt[static_cast<std::size_t>(sidx)].size());
    for (const int k : cache.sendAt[static_cast<std::size_t>(sidx)]) {
      const int frame = work.share.local[static_cast<std::size_t>(k)];
      const auto &p = cloud.pts[static_cast<std::size_t>(frame)];
      bag.push_back(GhostAtom{frame, p.atomID, p.type, p.molID, p.x, p.y, p.z});
    }
  }
  std::vector<int> sendCount(static_cast<std::size_t>(std::max(degree, 1)), 0);
  std::vector<int> sendDispl(sendCount.size(), 0);
  int sendN = 0;
  for (int sidx = 0; sidx < degree; sidx++) {
    sendCount[static_cast<std::size_t>(sidx)] =
        static_cast<int>(bags[static_cast<std::size_t>(sidx)].size());
    sendDispl[static_cast<std::size_t>(sidx)] = sendN;
    sendN += sendCount[static_cast<std::size_t>(sidx)];
  }
  std::vector<GhostAtom> send(static_cast<std::size_t>(std::max(sendN, 1)));
  for (int sidx = 0; sidx < degree; sidx++) {
    std::copy(bags[static_cast<std::size_t>(sidx)].begin(),
              bags[static_cast<std::size_t>(sidx)].end(),
              send.begin() + sendDispl[static_cast<std::size_t>(sidx)]);
  }
  const int atomBytes = static_cast<int>(sizeof(GhostAtom));
  std::vector<int> sendBytes(sendCount.size(), 0);
  std::vector<int> sendByteDispl(sendCount.size(), 0);
  for (int sidx = 0; sidx < degree; sidx++) {
    sendBytes[static_cast<std::size_t>(sidx)] =
        sendCount[static_cast<std::size_t>(sidx)] * atomBytes;
    sendByteDispl[static_cast<std::size_t>(sidx)] =
        sendDispl[static_cast<std::size_t>(sidx)] * atomBytes;
  }
  std::vector<int> recvBytes(sendCount.size(), 0);
  std::vector<int> recvByteDispl(sendCount.size(), 0);
  int sink = 0;
  const int *sendBytesPtr = degree == 0 ? &sink : sendBytes.data();
  const int *sendByteDisplPtr = degree == 0 ? &sink : sendByteDispl.data();
  MPI_Barrier(comm);
  const auto t0 = std::chrono::steady_clock::now();
  MPI_Neighbor_alltoall(sendBytesPtr, 1, MPI_INT, recvBytes.data(), 1, MPI_INT, graph);
  int recvN = 0;
  for (int sidx = 0; sidx < degree; sidx++) {
    recvByteDispl[static_cast<std::size_t>(sidx)] = recvN;
    recvN += recvBytes[static_cast<std::size_t>(sidx)];
  }
  std::vector<char> recv(static_cast<std::size_t>(std::max(recvN, 1)));
  MPI_Neighbor_alltoallv(send.data(), sendBytesPtr, sendByteDisplPtr, MPI_BYTE,
                         recv.data(), recvBytes.data(), recvByteDispl.data(), MPI_BYTE,
                         graph);
  MPI_Barrier(comm);
  work.commMs = std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - t0)
                    .count();
  (void)graph;

  // Place ghosts by frame index. The local list is sorted, so each
  // received atom is a binary search, and the point is rewritten in place.
  const bool fresh = work.cloud.pts.size() != work.share.local.size();
  if (fresh) {
    work.cloud.currentFrame = cloud.currentFrame;
    work.cloud.box = cloud.box;
    work.cloud.boxLow = cloud.boxLow;
    work.cloud.pts.assign(work.share.local.size(), {});
    work.cloud.idIndexMap.clear();
    work.cloud.idIndexMap.reserve(work.share.local.size());
  }
  gPlan.recvCount.resize(static_cast<std::size_t>(std::max(degree, 0)));
  int recvAtoms = 0;
  for (int sidx = 0; sidx < degree; sidx++) {
    gPlan.recvCount[static_cast<std::size_t>(sidx)] =
        recvBytes[static_cast<std::size_t>(sidx)] / std::max(atomBytes, 1);
    recvAtoms += gPlan.recvCount[static_cast<std::size_t>(sidx)];
  }
  const bool recordRecv = static_cast<int>(gPlan.recvLocal.size()) != recvAtoms;
  if (recordRecv) {
    gPlan.recvLocal.clear();
    gPlan.recvLocal.reserve(static_cast<std::size_t>(std::max(recvAtoms, 0)));
  }
  int misses = 0;
  for (std::size_t k = 0; k < work.share.local.size(); k++) {
    if (!work.share.owned[k]) {
      continue;
    }
    const auto &src = cloud.pts[static_cast<std::size_t>(work.share.local[k])];
    auto &dst = work.cloud.pts[k];
    dst.x = src.x;
    dst.y = src.y;
    dst.z = src.z;
    if (fresh) {
      dst.atomID = src.atomID;
      dst.molID = src.molID;
      dst.type = src.type;
      work.cloud.idIndexMap[dst.atomID] = static_cast<int>(k);
    }
  }
  for (int at = 0; at + atomBytes <= recvN; at += atomBytes) {
    GhostAtom atom;
    std::memcpy(&atom, recv.data() + at, sizeof(GhostAtom));
    const auto it = std::lower_bound(work.share.local.begin(), work.share.local.end(),
                                     atom.index);
    if (it == work.share.local.end() || *it != atom.index) {
      misses++;
      continue;
    }
    const std::size_t k = static_cast<std::size_t>(it - work.share.local.begin());
    auto &dst = work.cloud.pts[k];
    dst.x = atom.x;
    dst.y = atom.y;
    dst.z = atom.z;
    if (fresh) {
      dst.atomID = atom.atomID;
      dst.molID = atom.molID;
      dst.type = atom.type;
      work.cloud.idIndexMap[dst.atomID] = static_cast<int>(k);
    }
    if (recordRecv) {
      gPlan.recvLocal.push_back(static_cast<int>(k));
    }
  }
  if (fresh) {
    for (std::size_t k = 0; k < work.share.local.size(); k++) {
      if (work.share.owned[k] || work.cloud.pts[k].atomID != 0) {
        continue;
      }
      misses++;
      const auto &src = cloud.pts[static_cast<std::size_t>(work.share.local[k])];
      auto &dst = work.cloud.pts[k];
      dst = src;
      work.cloud.idIndexMap[dst.atomID] = static_cast<int>(k);
    }
  }
  work.cloud.nop = static_cast<int>(work.cloud.pts.size());
  int missed = 0;
  MPI_Allreduce(&misses, &missed, 1, MPI_INT, MPI_SUM, comm);
  if (rank == 0 && missed != 0) {
    std::fprintf(stderr, "ghost exchange missed %d atoms\n", missed);
  }
  return work;
}

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

void postBoundaryValues(double *values, int stride, MPI_Comm comm) {
  (void)comm;
  gValues = values;
  gStride = stride;
  gInflight = false;
  if (values == nullptr || stride <= 0 || gPlan.degree <= 0 ||
      gPlan.graph == MPI_COMM_NULL) {
    return;
  }
  const int nSend = static_cast<int>(gPlan.sendLocal.size());
  gSendPack.resize(static_cast<std::size_t>(std::max(nSend, 1) * stride));
  for (int i = 0; i < nSend; i++) {
    const int local = gPlan.sendLocal[static_cast<std::size_t>(i)];
    std::copy_n(values + static_cast<std::ptrdiff_t>(local * stride), stride,
                gSendPack.begin() + static_cast<std::ptrdiff_t>(i * stride));
  }
  gSendDoubles.resize(static_cast<std::size_t>(gPlan.degree));
  gSendDoubleDispl.resize(gSendDoubles.size());
  gRecvDoubles.resize(gSendDoubles.size());
  gRecvDoubleDispl.resize(gSendDoubles.size());
  int sendD = 0;
  int recvD = 0;
  for (int s = 0; s < gPlan.degree; s++) {
    gSendDoubles[static_cast<std::size_t>(s)] =
        gPlan.sendCount[static_cast<std::size_t>(s)] * stride;
    gSendDoubleDispl[static_cast<std::size_t>(s)] = sendD;
    sendD += gSendDoubles[static_cast<std::size_t>(s)];
    const int atoms = s < static_cast<int>(gPlan.recvCount.size())
                          ? gPlan.recvCount[static_cast<std::size_t>(s)]
                          : 0;
    gRecvDoubles[static_cast<std::size_t>(s)] = atoms * stride;
    gRecvDoubleDispl[static_cast<std::size_t>(s)] = recvD;
    recvD += gRecvDoubles[static_cast<std::size_t>(s)];
  }
  gRecvPack.resize(static_cast<std::size_t>(std::max(recvD, stride)));
  MPI_Ineighbor_alltoallv(gSendPack.data(), gSendDoubles.data(), gSendDoubleDispl.data(),
                          MPI_DOUBLE, gRecvPack.data(), gRecvDoubles.data(),
                          gRecvDoubleDispl.data(), MPI_DOUBLE, gPlan.graph, &gValueReq);
  gInflight = true;
}

void waitBoundaryValues() {
  if (!gInflight) {
    return;
  }
  const auto t0 = std::chrono::steady_clock::now();
  MPI_Wait(&gValueReq, MPI_STATUS_IGNORE);
  gValueCommMs += std::chrono::duration<double, std::milli>(
                      std::chrono::steady_clock::now() - t0)
                      .count();
  gInflight = false;
  const int stride = gStride;
  const int nRecv = static_cast<int>(gPlan.recvLocal.size());
  if (gValues == nullptr || stride <= 0) {
    return;
  }
  for (int i = 0; i < nRecv; i++) {
    std::copy_n(gRecvPack.begin() + static_cast<std::ptrdiff_t>(i * stride), stride,
                gValues + static_cast<std::ptrdiff_t>(gPlan.recvLocal[static_cast<std::size_t>(i)] *
                                                      stride));
  }
}

double takeValueCommMs() {
  const double ms = gValueCommMs;
  gValueCommMs = 0.0;
  return ms;
}

namespace {

int ownerOfFrame(const Cloud &cloud, int frame, const BlockGrid &grid) {
  const gen::FracBox b = gen::makeFracBox(cloud);
  const auto &p = cloud.pts[static_cast<std::size_t>(frame)];
  const double sz = (p.z - b.oz) / b.lz;
  const double sy = (p.y - b.oy - b.yz * sz) / b.ly;
  const double sx = (p.x - b.ox - b.xy * sy - b.xz * sz) / b.lx;
  const int ix = cellAlong(sx, grid.n[0]);
  const int iy = cellAlong(sy, grid.n[1]);
  const int iz = cellAlong(sz, grid.n[2]);
  return ix + grid.n[0] * (iy + grid.n[1] * iz);
}

} // namespace

std::vector<std::vector<int>> distributedRings(const Cloud &cloud, double cutoff,
                                               int maxDepth, MPI_Comm comm) {
  int rank = 0;
  int nRanks = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nRanks);
  if (nRanks <= 1 || cloud.nop == 0) {
    std::vector<char> all(static_cast<std::size_t>(std::max(cloud.nop, 0)), 1);
    return ringsOf(cloud, cutoff, maxDepth, all);
  }
  const gen::FracBox cell = gen::makeFracBox(cloud);
  if (std::fabs(cell.xy) + std::fabs(cell.xz) + std::fabs(cell.yz) > 1.0e-8) {
    // A skewed shell does not fit in the block-distance halo. The rank
    // already holds the frame, so it reads those rows directly.
    return rings(cloud, cutoff, maxDepth, rank, nRanks);
  }
  // ringHaloBonds(maxDepth) is floor(maxDepth/2). The coordinate skin is
  // that many cutoffs: every primitive ring through an owned atom, and
  // every shortcut of it, stays inside the skin.
  const Work &work = gatherWork(cloud, ringHaloBonds(maxDepth) * cutoff, comm);
  auto local = ownedRings(work.cloud, cutoff, maxDepth, work.share.owned);
  for (auto &ring : local) {
    for (int &v : ring) {
      v = work.share.local[static_cast<std::size_t>(v)];
    }
  }
  return local;
}


RingSummary reduceRings(const Cloud &cloud, double cutoff, int maxDepth,
                        MPI_Comm comm, int root) {
  int rank = 0;
  int nRanks = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nRanks);
  const int nAtoms = static_cast<int>(cloud.pts.size());
  const auto mine = distributedRings(cloud, cutoff, maxDepth, comm);

  RingSummary s;
  std::vector<long long> count(static_cast<std::size_t>(std::max(maxDepth, 0)) + 1, 0);
  for (const auto &ring : mine) {
    if (ring.size() < count.size()) {
      count[ring.size()]++;
    }
  }
  s.sizeCount.assign(count.size(), 0);
  MPI_Allreduce(count.data(), s.sizeCount.data(), static_cast<int>(count.size()),
                MPI_LONG_LONG, MPI_SUM, comm);

  std::vector<int> type(static_cast<std::size_t>(std::max(nAtoms, 0)), 1);
  for (const auto &ring : mine) {
    const int size = static_cast<int>(ring.size());
    for (const int v : ring) {
      if (v >= 0 && v < nAtoms) {
        type[static_cast<std::size_t>(v)] =
            std::max(type[static_cast<std::size_t>(v)], size);
      }
    }
  }
  const gen::FracBox box = gen::makeFracBox(cloud);
  const double w[3] = {box.wx, box.wy, box.wz};
  const BlockGrid grid = chooseGrid(nRanks, w, cutoff);
  std::vector<int> owner(static_cast<std::size_t>(std::max(nAtoms, 0)), 0);
  if (nRanks > 1 && box.ok) {
    for (int i = 0; i < nAtoms; i++) {
      owner[static_cast<std::size_t>(i)] = ownerOfFrame(cloud, i, grid);
    }
  }
  std::vector<int> sendCount(static_cast<std::size_t>(std::max(nRanks, 1)), 0);
  for (int i = 0; i < nAtoms; i++) {
    const int o = owner[static_cast<std::size_t>(i)];
    if (type[static_cast<std::size_t>(i)] > 1 && o != rank && o >= 0 && o < nRanks) {
      sendCount[static_cast<std::size_t>(o)] += 2;
    }
  }
  std::vector<int> sendDispl(static_cast<std::size_t>(nRanks), 0);
  std::exclusive_scan(sendCount.begin(), sendCount.end(), sendDispl.begin(), 0);
  const int sendN = nRanks == 0 ? 0 : sendDispl.back() + sendCount.back();
  std::vector<int> send(static_cast<std::size_t>(std::max(sendN, 1)));
  {
    std::vector<int> at = sendDispl;
    for (int i = 0; i < nAtoms; i++) {
      const int o = owner[static_cast<std::size_t>(i)];
      if (type[static_cast<std::size_t>(i)] > 1 && o != rank && o >= 0 && o < nRanks) {
        int &pos = at[static_cast<std::size_t>(o)];
        send[static_cast<std::size_t>(pos++)] = i;
        send[static_cast<std::size_t>(pos++)] = type[static_cast<std::size_t>(i)];
      }
    }
  }
  std::vector<int> recvCount(static_cast<std::size_t>(std::max(nRanks, 1)), 0);
  MPI_Alltoall(sendCount.data(), 1, MPI_INT, recvCount.data(), 1, MPI_INT, comm);
  std::vector<int> recvDispl(static_cast<std::size_t>(nRanks), 0);
  std::exclusive_scan(recvCount.begin(), recvCount.end(), recvDispl.begin(), 0);
  const int recvN = nRanks == 0 ? 0 : recvDispl.back() + recvCount.back();
  std::vector<int> recv(static_cast<std::size_t>(std::max(recvN, 1)));
  MPI_Alltoallv(send.data(), sendCount.data(), sendDispl.data(), MPI_INT, recv.data(),
                recvCount.data(), recvDispl.data(), MPI_INT, comm);
  for (int j = 0; j + 1 < recvN; j += 2) {
    const int frame = recv[static_cast<std::size_t>(j)];
    if (frame >= 0 && frame < nAtoms) {
      type[static_cast<std::size_t>(frame)] =
          std::max(type[static_cast<std::size_t>(frame)], recv[static_cast<std::size_t>(j + 1)]);
    }
  }

  std::vector<int> owned;
  for (int i = 0; i < nAtoms; i++) {
    if (owner[static_cast<std::size_t>(i)] == rank && type[static_cast<std::size_t>(i)] > 1) {
      owned.push_back(i);
      owned.push_back(type[static_cast<std::size_t>(i)]);
    }
  }
  const int ownedCount = static_cast<int>(owned.size());
  std::vector<int> counts(rank == root ? static_cast<std::size_t>(nRanks) : 0);
  MPI_Gather(&ownedCount, 1, MPI_INT, counts.data(), 1, MPI_INT, root, comm);
  std::vector<int> displs(counts.size(), 0);
  std::vector<int> all;
  if (rank == root) {
    std::exclusive_scan(counts.begin(), counts.end(), displs.begin(), 0);
    all.resize(static_cast<std::size_t>(displs.back() + counts.back()));
  }
  MPI_Gatherv(owned.data(), ownedCount, MPI_INT, all.data(), counts.data(),
              displs.data(), MPI_INT, root, comm);
  if (rank == root) {
    s.atomType.assign(static_cast<std::size_t>(nAtoms), 1);
    for (std::size_t j = 0; j + 1 < all.size(); j += 2) {
      s.atomType[static_cast<std::size_t>(all[j])] = all[j + 1];
    }
  }
  return s;
}
#endif

} // namespace seams::domain
