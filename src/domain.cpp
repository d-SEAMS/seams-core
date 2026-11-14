#include <domain.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <unordered_map>

#ifdef SEAMS_HAS_LINKCELL
#include <linkcell.h>
// A library from before lc_owned_mask still links. The header's copy wins.
extern "C" __attribute__((weak)) int lc_owned_mask(const int *nxyz, const int *assign, std::size_t ncell,
                             const int *point_cell, std::size_t nop, int rank,
                             const std::uint32_t *hwidth, int *hit) {
  (void)nxyz;
  (void)assign;
  (void)ncell;
  (void)point_cell;
  (void)nop;
  (void)rank;
  (void)hwidth;
  (void)hit;
  return 1;
}
#endif

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

#ifdef SEAMS_HAS_LINKCELL
// Point-count cut of the periodic cell grid. Measured weights replace it
// on the next call of the same analysis. The grid itself is reused.
struct CountCut {
  int nop = -1;
  int parts = -1;
  double halo = -1.0;
  double lx = 0, ly = 0, lz = 0;
  double origin[3] = {0, 0, 0};
  int nxyz[3] = {0, 0, 0};
  double width[3] = {0, 0, 0};
  std::uint32_t hwidth[3] = {0, 0, 0};
  std::vector<int> assign;
  std::vector<int> pointCell;
};

CountCut gCount;
bool gForceBlock = false;
// 0 chooses from the measured block loads, 1 is the block cut, 2 is the cell cut.
int gCutMode = 0;

struct MeasuredCut {
  int nop = -1;
  int parts = -1;
  double halo = -1.0;
  bool ready = false;
  std::vector<std::uint64_t> weight;
  std::vector<int> assign;
};

MeasuredCut gMeasured[5];

bool sameBox(const CountCut &cut, const Cloud &cloud, double halo, int parts) {
  return cut.nop == cloud.nop && cut.parts == parts && cut.halo == halo &&
         cloud.box.size() == 3 && cut.lx == cloud.box[0] && cut.ly == cloud.box[1] &&
         cut.lz == cloud.box[2];
}

bool buildCountCut(const Cloud &cloud, double halo, int parts) {
  if (sameBox(gCount, cloud, halo, parts) &&
      gCount.pointCell.size() == cloud.pts.size() && !gCount.assign.empty()) {
    return true;
  }
  if (cloud.box.size() != 3 || !(halo > 0.0) || parts < 1 || cloud.nop <= 0) {
    return false;
  }
  std::vector<double> xyz(static_cast<std::size_t>(cloud.nop) * 3);
  for (int i = 0; i < cloud.nop; i++) {
    const auto &p = cloud.pts[static_cast<std::size_t>(i)];
    xyz[static_cast<std::size_t>(i) * 3] = p.x;
    xyz[static_cast<std::size_t>(i) * 3 + 1] = p.y;
    xyz[static_cast<std::size_t>(i) * 3 + 2] = p.z;
  }
  lc_cell box = lc_cell_ortho(cloud.box[0], cloud.box[1], cloud.box[2]);
  if (cloud.boxLow.size() >= 3) {
    box.ox = cloud.boxLow[0];
    box.oy = cloud.boxLow[1];
    box.oz = cloud.boxLow[2];
  }
  std::size_t ncell = 0;
  // Bins finer than the cutoff let the distance filter drop the outer
  // part of a painted cell instead of keeping the whole neighbour bin.
  const double hint = halo * 0.5;
  const int query = lc_mesh_cut(xyz.data(), static_cast<std::size_t>(cloud.nop), &box, halo,
                                 nullptr, hint, static_cast<std::size_t>(parts), 2, 0.0, 0, nullptr,
                                 nullptr, nullptr, nullptr, 0, 0, &ncell);
  if (query != 0 || ncell == 0) {
    return false;
  }
  std::vector<int> assign(ncell, 0);
  std::vector<int> point(static_cast<std::size_t>(cloud.nop), 0);
  int nxyz[3] = {0, 0, 0};
  double width[3] = {0, 0, 0};
  const int filled =
      lc_mesh_cut(xyz.data(), static_cast<std::size_t>(cloud.nop), &box, halo, nullptr, hint,
                  static_cast<std::size_t>(parts), 2, 0.0, 0, nxyz, width, assign.data(),
                  point.data(), ncell, static_cast<std::size_t>(cloud.nop), &ncell);
  if (filled != 0) {
    return false;
  }
  std::uint32_t hwidth[3] = {0, 0, 0};
  if (lc_halo_width(halo, 1, width, hwidth) != 0) {
    return false;
  }
  gCount.nop = cloud.nop;
  gCount.parts = parts;
  gCount.halo = halo;
  gCount.lx = cloud.box[0];
  gCount.ly = cloud.box[1];
  gCount.lz = cloud.box[2];
  gCount.origin[0] = cloud.boxLow.size() > 0 ? cloud.boxLow[0] : 0.0;
  gCount.origin[1] = cloud.boxLow.size() > 1 ? cloud.boxLow[1] : 0.0;
  gCount.origin[2] = cloud.boxLow.size() > 2 ? cloud.boxLow[2] : 0.0;
  for (int d = 0; d < 3; d++) {
    gCount.nxyz[d] = nxyz[d];
    gCount.width[d] = width[d];
    gCount.hwidth[d] = hwidth[d];
  }
  gCount.assign = std::move(assign);
  gCount.pointCell = std::move(point);
  return true;
}

bool shareFromAssign(const Cloud &cloud, int rank, const std::vector<int> &assign, Share &share) {
  if (gCount.pointCell.size() != cloud.pts.size() || assign.size() != gCount.assign.size()) {
    return false;
  }
  const int ncell = static_cast<int>(assign.size());
  std::vector<int> hit(static_cast<std::size_t>(cloud.nop), 0);
  // rank_mask paints the stencil while it walks the cells this rank owns.
  if (lc_owned_mask(gCount.nxyz, assign.data(), assign.size(), gCount.pointCell.data(),
                    static_cast<std::size_t>(cloud.nop), rank, gCount.hwidth, hit.data()) != 0) {
    return false;
  }
  share.local.clear();
  share.owned.clear();
  share.owner.clear();
  for (int i = 0; i < cloud.nop; i++) {
    const int flag = hit[static_cast<std::size_t>(i)];
    if (flag == 0) {
      continue;
    }
    const int cell = gCount.pointCell[static_cast<std::size_t>(i)];
    share.local.push_back(i);
    share.owned.push_back(flag == 2 ? 1 : 0);
    share.owner.push_back(cell >= 0 && cell < ncell ? assign[static_cast<std::size_t>(cell)] : rank);
  }
  return true;
}

void linkcellSends(const Cloud &cloud, const Share &share, const std::vector<int> &assign,
                   const std::vector<int> &nbrs, std::vector<std::vector<int>> &sendAt) {
  sendAt.assign(nbrs.size(), {});
  std::vector<int> hit(static_cast<std::size_t>(cloud.nop), 0);
  for (std::size_t sidx = 0; sidx < nbrs.size(); sidx++) {
    const int other = nbrs[sidx];
    if (lc_owned_mask(gCount.nxyz, assign.data(), assign.size(), gCount.pointCell.data(),
                      static_cast<std::size_t>(cloud.nop), other, gCount.hwidth, hit.data()) != 0) {
      continue;
    }
    auto &dest = sendAt[sidx];
    for (std::size_t k = 0; k < share.local.size(); k++) {
      if (!share.owned[k]) {
        continue;
      }
      const int frame = share.local[k];
      if (frame >= 0 && static_cast<std::size_t>(frame) < hit.size() &&
          hit[static_cast<std::size_t>(frame)] == 1) {
        dest.push_back(static_cast<int>(k));
      }
    }
  }
}
#endif

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

#ifdef SEAMS_HAS_LINKCELL
  if (!gForceBlock && cloud.box.size() == 3 && buildCountCut(cloud, halo, nRanks)) {
    Share cut;
    if (shareFromAssign(cloud, rank, gCount.assign, cut)) {
      return cut;
    }
  }
#endif
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
    const auto &p = cloud.pts[static_cast<std::size_t>(share.local[k])];
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
seams::domain::RingPhases gRingPhases;


//! Bit mix of the coordinates. A later frame with the same atom ids still
//! rebuilds its rows when a point has moved.
std::uint64_t coordinateKey(const Cloud &local) {
  std::uint64_t h = 1469598103934665603ull;
  for (int i = 0; i < local.nop; i++) {
    const auto &p = local.pts[static_cast<std::size_t>(i)];
    const double xyz[3] = {p.x, p.y, p.z};
    for (const double s : xyz) {
      std::uint64_t u = 0;
      std::memcpy(&u, &s, sizeof u);
      h ^= u;
      h *= 0x100000001b3ull;
    }
  }
  return h;
}

//! Rings on a cloud already limited to one share, by local index.
const std::vector<std::vector<int>> &ringsOf(const Cloud &local, double cutoff,
                                                  int maxDepth, const std::vector<char> &owned) {
  if (local.nop == 0) {
    static const std::vector<std::vector<int>> empty;
    return empty;
  }
  // Orthorhombic frames take the cutoff CSR directly. The ring search
  // reads those offsets and columns; it does not expand them into a
  // vector of rows first.
  {
    std::vector<int> offsets;
    std::vector<int> cols;
    if (nneigh::orthoFullCsr(local, nullptr, cutoff, offsets, cols)) {
      thread_local std::vector<std::vector<int>> csrRings;
      const auto search0 = std::chrono::steady_clock::now();
      csrRings = primitive::ringNetwork(offsets, cols, maxDepth, owned);
      const double searchAll =
          std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - search0)
              .count();
      gRingPhases.primitiveMs = primitive::takePrimitiveMs();
      gRingPhases.searchMs = std::max(0.0, searchAll - gRingPhases.primitiveMs);
      gRingPhases.rebuildMs = 0;
      gRingPhases.reduceMs = 0;
      return csrRings;
    }
  }
  // Rows are reused only when the points themselves are unchanged. Two
  // slots cover the block and cell clouds the bench alternates. A new
  // frame keeps the domain cut and rebuilds these rows.
  struct RowCache {
    int nop = -1;
    double cutoff = 0;
    int maxDepth = -1;
    unsigned stamp = 0;
    std::uint64_t coord = 0;
    std::vector<int> id;
    std::vector<char> owned;
    std::vector<std::vector<int>> raw;
    std::vector<std::vector<int>> trimmed;
    std::vector<std::vector<int>> rings;
  };
  // One cache per calling thread. Overlapping rings() calls used to
  // share these rows and tear each other's neighbour list.
  thread_local RowCache cache[2];
  thread_local unsigned clock = 1;
  const std::uint64_t coord = coordinateKey(local);
  int hit = -1;
  for (int s = 0; s < 2; s++) {
    const RowCache &c = cache[s];
    if (c.nop != local.nop || c.cutoff != cutoff || c.coord != coord ||
        static_cast<int>(c.id.size()) != local.nop) {
      continue;
    }
    bool same = true;
    for (int i = 0; i < local.nop; i++) {
      if (c.id[static_cast<std::size_t>(i)] != local.pts[static_cast<std::size_t>(i)].atomID) {
        same = false;
        break;
      }
    }
    if (same) {
      hit = s;
      break;
    }
  }
  if (hit >= 0 && cache[hit].maxDepth != -1 && cache[hit].maxDepth != maxDepth) {
    cache[hit].nop = -1;
    hit = -1;
  }
  double rebuildMs = 0;
  if (hit < 0) {
    const auto rebuild0 = std::chrono::steady_clock::now();
    if (cache[0].nop < 0) {
      hit = 0;
    } else if (cache[1].nop < 0) {
      hit = 1;
    } else {
      hit = cache[0].stamp <= cache[1].stamp ? 0 : 1;
    }
    RowCache &c = cache[hit];
    // A partial share is the linkcell cut. Fold the window onto the owned
    // atoms and give every atom in that window its partners. One rank
    // owns the frame, so it keeps the full-box list.
    int nMarked = 0;
    for (const char bit : owned) {
      nMarked += bit != 0;
    }
    const bool partial = !owned.empty() && nMarked > 0 && nMarked < local.nop;
    bool built = false;
    if (partial) {
      built = nneigh::blockWindowRows(local, cutoff, &owned, c.raw, false, true);
    }
    if (!built && !nneigh::blockWindowRows(local, cutoff, nullptr, c.raw, false)) {
      if (!nneigh::refillIndexRows(local, cutoff, c.raw)) {
        c.raw = nneigh::getNewNeighbourListByIndex(local, cutoff);
      }
      built = false;
    }
    if (built) {
      const int nBuilt = static_cast<int>(c.raw.size());
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nBuilt >= kParallelMin && !omp_in_parallel())
#endif
      for (int k = 0; k < nBuilt; k++) {
        auto &row = c.raw[static_cast<std::size_t>(k)];
        if (row.size() < 3) {
          continue;
        }
        std::size_t w = 2;
        for (std::size_t j = 2; j < row.size(); j++) {
          if (row[j] != row[w - 1]) {
            row[w++] = row[j];
          }
        }
        row.resize(w);
      }
    }
    if (!built) {
      const int nBuilt = static_cast<int>(c.raw.size());
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nBuilt >= kParallelMin && !omp_in_parallel())
#endif
      for (int k = 0; k < nBuilt; k++) {
        if (c.raw[static_cast<std::size_t>(k)].size() > 2) {
          std::sort(c.raw[static_cast<std::size_t>(k)].begin() + 1,
                    c.raw[static_cast<std::size_t>(k)].end());
        }
      }
    }
    c.nop = local.nop;
    c.cutoff = cutoff;
    c.coord = coord;
    c.maxDepth = -1;
    c.id.resize(static_cast<std::size_t>(local.nop));
    for (int i = 0; i < local.nop; i++) {
      c.id[static_cast<std::size_t>(i)] = local.pts[static_cast<std::size_t>(i)].atomID;
    }
    rebuildMs += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() -
                                                           rebuild0)
                     .count();
  }
  cache[hit].stamp = clock++;
  const bool sameTrim = cache[hit].maxDepth == maxDepth && cache[hit].owned == owned;
  const auto trim0 = std::chrono::steady_clock::now();
  auto &nList = cache[hit].raw;
  if (!sameTrim) {
    cache[hit].owned = owned;
    cache[hit].maxDepth = maxDepth;
  }
  const int nLocal = static_cast<int>(nList.size());
  // Halo.lean: a primitive ring of length at most k lies in the closed
  // bond ball of radius floor(k/2) about an owned vertex. Vertices past
  // that ball are not members and are not shortcuts.
  if (!sameTrim && !owned.empty() && nLocal > 0) {
    const int radius = ringHops(maxDepth);
    std::vector<int> dist(static_cast<std::size_t>(nLocal), -1);
    std::vector<int> queue;
    queue.reserve(static_cast<std::size_t>(nLocal));
    const int nOwned = std::min(nLocal, static_cast<int>(owned.size()));
    for (int i = 0; i < nOwned; i++) {
      if (!owned[static_cast<std::size_t>(i)]) {
        continue;
      }
      dist[static_cast<std::size_t>(i)] = 0;
      queue.push_back(i);
    }
    for (std::size_t head = 0; head < queue.size(); head++) {
      const int u = queue[head];
      if (dist[static_cast<std::size_t>(u)] >= radius) {
        continue;
      }
      const auto &row = nList[static_cast<std::size_t>(u)];
      for (std::size_t j = 1; j < row.size(); j++) {
        const int v = row[j];
        if (v < 0 || v >= nLocal || dist[static_cast<std::size_t>(v)] >= 0) {
          continue;
        }
        dist[static_cast<std::size_t>(v)] = dist[static_cast<std::size_t>(u)] + 1;
        queue.push_back(v);
      }
    }
    for (int i = 0; i < nLocal; i++) {
      auto &row = nList[static_cast<std::size_t>(i)];
      if (dist[static_cast<std::size_t>(i)] < 0) {
        row.resize(std::min<std::size_t>(row.size(), 1));
        continue;
      }
      std::size_t w = 1;
      for (std::size_t j = 1; j < row.size(); j++) {
        const int v = row[j];
        if (v >= 0 && v < nLocal && dist[static_cast<std::size_t>(v)] >= 0) {
          row[w++] = v;
        }
      }
      row.resize(w);
    }
  }
  if (!sameTrim) {
    rebuildMs += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - trim0)
                     .count();
  }
  const std::vector<std::vector<int>> &graph = nList;
  const std::vector<char> &sources = cache[hit].owned.empty() ? owned : cache[hit].owned;
  const auto search0 = std::chrono::steady_clock::now();
  auto found = primitive::ringNetwork(graph, maxDepth, sources);
  const double searchAll =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - search0).count();
  gRingPhases.primitiveMs = primitive::takePrimitiveMs();
  gRingPhases.searchMs = std::max(0.0, searchAll - gRingPhases.primitiveMs);
  gRingPhases.rebuildMs = rebuildMs;
  gRingPhases.reduceMs = 0;
  cache[hit].rings = std::move(found);
  return cache[hit].rings;
}

} // namespace

RingPhases lastRingPhases() { return gRingPhases; }

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

const std::vector<std::vector<int>> &ownedRings(const Cloud &local, double cutoff,
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

void setCutMode(int mode) {
#ifdef SEAMS_HAS_LINKCELL
  gCutMode = mode;
#else
  (void)mode;
#endif
}

namespace {
void (*gHaloHook)(void *) = nullptr;
void *gHaloHookCtx = nullptr;
} // namespace

void setHaloHook(void (*fn)(void *), void *ctx) {
  gHaloHook = fn;
  gHaloHookCtx = ctx;
}

void markInterior(const Work &work, double cutoff, std::vector<char> &interior) {
  interior.assign(work.cloud.pts.size(), 0);
  const gen::FracBox box = gen::makeFracBox(work.cloud);
  if (!box.ok || !(cutoff > 0.0) || work.share.owned.empty()) {
    return;
  }
  const double len[3] = {box.lx, box.ly, box.lz};
  double axis[3][2] = {{1.0e300, -1.0e300}, {1.0e300, -1.0e300}, {1.0e300, -1.0e300}};
  std::vector<double> coord;
  coord.reserve(work.share.owned.size());
  for (int d = 0; d < 3; d++) {
    coord.clear();
    for (std::size_t k = 0; k < work.share.owned.size(); k++) {
      if (!work.share.owned[k] || k >= work.cloud.pts.size()) {
        continue;
      }
      const auto &p = work.cloud.pts[k];
      const double c[3] = {p.x, p.y, p.z};
      coord.push_back(c[d]);
    }
    if (coord.empty()) {
      return;
    }
    std::sort(coord.begin(), coord.end());
    const double length = len[d];
    double bestGap = coord.front() + length - coord.back();
    int after = static_cast<int>(coord.size()) - 1;
    for (int i = 0; i + 1 < static_cast<int>(coord.size()); i++) {
      const double gap = coord[static_cast<std::size_t>(i + 1)] - coord[static_cast<std::size_t>(i)];
      if (gap > bestGap) {
        bestGap = gap;
        after = i;
      }
    }
    const double start =
        after + 1 == static_cast<int>(coord.size()) ? coord.front() : coord[static_cast<std::size_t>(after + 1)];
    axis[d][0] = start;
    axis[d][1] = start + std::max(length - bestGap, 0.0);
  }
  for (std::size_t k = 0; k < work.share.owned.size() && k < work.cloud.pts.size(); k++) {
    if (!work.share.owned[k]) {
      continue;
    }
    const auto &p = work.cloud.pts[k];
    const double c[3] = {p.x, p.y, p.z};
    bool inside = true;
    for (int d = 0; d < 3 && inside; d++) {
      double s = c[d];
      if (s < axis[d][0]) {
        s += len[d];
      }
      if (s > axis[d][1]) {
        s -= len[d];
      }
      if (!(s > axis[d][0] + cutoff && s < axis[d][1] - cutoff)) {
        inside = false;
      }
    }
    if (inside) {
      interior[k] = 1;
    }
  }
}

namespace {
#ifdef SEAMS_HAS_LINKCELL
// Heaviest block rank divided by the mean. The caller keeps the block
// cut while this stays under 1.025.
double blockLoadRatio(const std::vector<std::uint64_t> &weight, int nRanks, double halo) {
  if (weight.empty() || nRanks <= 1 || gCount.nxyz[0] <= 0) {
    return 1.0;
  }
  const double wlen[3] = {gCount.lx, gCount.ly, gCount.lz};
  const BlockGrid grid = chooseGrid(nRanks, wlen, halo);
  std::vector<double> load(static_cast<std::size_t>(nRanks), 0.0);
  const int nx = gCount.nxyz[0];
  const int ny = gCount.nxyz[1];
  const int nz = gCount.nxyz[2];
  double total = 0.0;
  for (std::size_t c = 0; c < weight.size(); c++) {
    if (weight[c] == 0) {
      continue;
    }
    const int ix = static_cast<int>(c) % nx;
    const int iy = (static_cast<int>(c) / nx) % ny;
    const int iz = static_cast<int>(c) / (nx * ny);
    const double s[3] = {(ix + 0.5) / nx, (iy + 0.5) / ny, (iz + 0.5) / nz};
    const int blk = cellAlong(s[0], grid.n[0]) +
                    grid.n[0] * (cellAlong(s[1], grid.n[1]) + grid.n[1] * cellAlong(s[2], grid.n[2]));
    if (blk >= 0 && blk < nRanks) {
      load[static_cast<std::size_t>(blk)] += static_cast<double>(weight[c]);
    }
    total += static_cast<double>(weight[c]);
  }
  if (!(total > 0.0)) {
    return 1.0;
  }
  double most = 0.0;
  for (const double v : load) {
    most = std::max(most, v);
  }
  return most / (total / nRanks);
}

double assignmentLoad(const std::vector<std::uint64_t> &weight, const std::vector<int> &assign,
                      int nRanks) {
  if (weight.size() != assign.size() || nRanks <= 1) {
    return 1.0;
  }
  std::vector<double> load(static_cast<std::size_t>(nRanks), 0.0);
  double total = 0.0;
  for (std::size_t c = 0; c < weight.size(); c++) {
    const int owner = assign[c];
    if (owner < 0 || owner >= nRanks || weight[c] == 0) {
      continue;
    }
    load[static_cast<std::size_t>(owner)] += static_cast<double>(weight[c]);
    total += static_cast<double>(weight[c]);
  }
  if (!(total > 0.0)) {
    return 1.0;
  }
  double most = 0.0;
  for (const double v : load) {
    most = std::max(most, v);
  }
  return most / (total / nRanks);
}

// ghosts/local on the rank with the thickest skin. The same count the
// bench reports as the halo fraction.
double blockGhostFraction(const Cloud &cloud, int nRanks, double halo) {
  const gen::FracBox box = gen::makeFracBox(cloud);
  if (!box.ok || nRanks <= 1 || !(halo > 0.0)) {
    return 0.0;
  }
  const double wlen[3] = {box.wx, box.wy, box.wz};
  const BlockGrid grid = chooseGrid(nRanks, wlen, halo);
  const double halo2 = halo * halo;
  std::vector<int> blocks(static_cast<std::size_t>(nRanks) * 3, 0);
  for (int rank = 0; rank < nRanks; rank++) {
    blockAt(rank, grid.n, blocks.data() + static_cast<std::ptrdiff_t>(rank * 3));
  }
  std::vector<int> owned(static_cast<std::size_t>(nRanks), 0);
  std::vector<int> ghost(static_cast<std::size_t>(nRanks), 0);
  for (int i = 0; i < cloud.nop; i++) {
    const auto &p = cloud.pts[static_cast<std::size_t>(i)];
    const double sz = (p.z - box.oz) / box.lz;
    const double sy = (p.y - box.oy - box.yz * sz) / box.ly;
    const double sx = (p.x - box.ox - box.xy * sy - box.xz * sz) / box.lx;
    const double s[3] = {sx, sy, sz};
    const int owner = cellAlong(sx, grid.n[0]) +
                      grid.n[0] * (cellAlong(sy, grid.n[1]) + grid.n[1] * cellAlong(sz, grid.n[2]));
    if (owner < 0 || owner >= nRanks) {
      continue;
    }
    owned[static_cast<std::size_t>(owner)]++;
    for (int rank = 0; rank < nRanks; rank++) {
      if (rank == owner) {
        continue;
      }
      if (blockDist2(s, blocks.data() + static_cast<std::ptrdiff_t>(rank * 3), grid.n, wlen) <=
          halo2) {
        ghost[static_cast<std::size_t>(rank)]++;
      }
    }
  }
  double worst = 0.0;
  for (int rank = 0; rank < nRanks; rank++) {
    const int local = owned[static_cast<std::size_t>(rank)] + ghost[static_cast<std::size_t>(rank)];
    if (local > 0) {
      worst = std::max(worst, static_cast<double>(ghost[static_cast<std::size_t>(rank)]) / local);
    }
  }
  return worst;
}

double cellGhostFraction(const Cloud &cloud, const std::vector<int> &assign, int nRanks,
                         double halo) {
  if (assign.size() != gCount.assign.size() || cloud.pts.size() != gCount.pointCell.size() ||
      nRanks <= 1) {
    return 0.0;
  }
  std::vector<int> hit(static_cast<std::size_t>(cloud.nop), 0);
  double worst = 0.0;
  for (int rank = 0; rank < nRanks; rank++) {
    if (lc_owned_mask(gCount.nxyz, assign.data(), assign.size(), gCount.pointCell.data(),
                      static_cast<std::size_t>(cloud.nop), rank, gCount.hwidth, hit.data()) != 0) {
      return 0.0;
    }
    int owned = 0;
    int ghost = 0;
    for (const int flag : hit) {
      if (flag == 2) {
        owned++;
      } else if (flag == 1) {
        ghost++;
      }
    }
    const int local = owned + ghost;
    if (local > 0) {
      worst = std::max(worst, static_cast<double>(ghost) / local);
    }
  }
  return worst;
}
#endif
} // namespace

Work &gatherWork(const Cloud &cloud, double halo, MPI_Comm comm, int workTag) {
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
    int tag = -1;
    int epoch = -1;
    double halo = -1.0;
    double haloFraction = 0.0;
    Share share;
    std::vector<int> nbrs;
    std::vector<std::vector<int>> sendAt;
    MPI_Comm graph = MPI_COMM_NULL;
  };
  static GraphCache byMode[2];
#ifdef SEAMS_HAS_LINKCELL
  if (const char *forced = std::getenv("SEAMS_FORCE_CUT")) {
    if (std::strcmp(forced, "cell") == 0) {
      gCutMode = 2;
    } else if (std::strcmp(forced, "block") == 0) {
      gCutMode = 1;
    }
  }
  const bool weighted = workTag > 0 && workTag < 5 && gMeasured[workTag].ready &&
                        gMeasured[workTag].nop == cloud.nop && gMeasured[workTag].parts == nRanks &&
                        gMeasured[workTag].halo == halo;
  const int epoch = weighted ? 1 : 0;
  double ratio = 1.0;
  double blockCost = 1.0;
  double cellCost = 1.0;
  bool autoCell = false;
  // An orthorhombic frame runs on the linkcell bisection. The geometric
  // blocks remain the fallback when that mesh cannot be built, and when a
  // caller forces them.
  bool useCell = gCutMode != 1 && cloud.box.size() == 3 && buildCountCut(cloud, halo, nRanks);
  if (useCell && std::getenv("SEAMS_CUT_LOG") != nullptr) {
    struct CostKey {
      int nop = -1;
      int parts = -1;
      int tag = -1;
      int epoch = -1;
      double halo = -1.0;
      bool autoCell = false;
      double ratio = 1.0;
      double blockCost = 1.0;
      double cellCost = 1.0;
      double blockGhost = 0.0;
      double cellGhost = 0.0;
    };
    static CostKey cached;
    if (cached.nop == cloud.nop && cached.parts == nRanks && cached.tag == workTag &&
        cached.epoch == epoch && cached.halo == halo) {
      ratio = cached.ratio;
      blockCost = cached.blockCost;
      cellCost = cached.cellCost;
      autoCell = cached.autoCell;
    } else {
      const std::vector<std::uint64_t> *weight = nullptr;
      std::vector<std::uint64_t> counts;
      const std::vector<int> *cellAssign = &gCount.assign;
      if (weighted) {
        weight = &gMeasured[workTag].weight;
        if (gMeasured[workTag].assign.size() == gCount.assign.size()) {
          cellAssign = &gMeasured[workTag].assign;
        }
      } else {
        counts.assign(gCount.assign.size(), 0);
        for (const int cell : gCount.pointCell) {
          if (cell >= 0 && static_cast<std::size_t>(cell) < counts.size()) {
            counts[static_cast<std::size_t>(cell)]++;
          }
        }
        weight = &counts;
      }
      ratio = blockLoadRatio(*weight, nRanks, halo);
      const double blockGhost = blockGhostFraction(cloud, nRanks, halo);
      const double cellLoad = assignmentLoad(*weight, *cellAssign, nRanks);
      const double cellGhost = cellGhostFraction(cloud, *cellAssign, nRanks, halo);
      // Load is the heaviest rank over the mean. Ghost fraction is the
      // thickest skin. Their sum is the cost the bench's wall follows:
      // a balanced block cut with a thick skin loses to a thin cell skin.
      blockCost = ratio + blockGhost;
      cellCost = cellLoad + cellGhost;
      autoCell = cellCost < blockCost;
      cached = CostKey{cloud.nop, nRanks, workTag, epoch, halo, autoCell, ratio, blockCost,
                       cellCost, blockGhost, cellGhost};
    }
  }
  if (rank == 0 && std::getenv("SEAMS_CUT_LOG") != nullptr) {
    std::fprintf(stderr, "cut %s ratio %.3f rule %s costb %.3f costc %.3f\n",
                 useCell ? "cell" : "block", ratio, autoCell ? "cell" : "block", blockCost,
                 cellCost);
  }
#else
  const bool weighted = false;
  const int epoch = 0;
  const bool useCell = false;
#endif
  GraphCache &cache = byMode[useCell ? 1 : 0];
  const bool reuse = cache.nop == cloud.nop && cache.nRanks == nRanks && cache.rank == rank &&
                     cache.halo == halo && cache.tag == workTag && cache.epoch == epoch;
  bool usedCut = false;
  if (!reuse) {
#ifdef SEAMS_HAS_LINKCELL
    if (useCell) {
      const std::vector<int> &assign = weighted ? gMeasured[workTag].assign : gCount.assign;
      usedCut = shareFromAssign(cloud, rank, assign, work.share);
    }
    if (!usedCut) {
      gForceBlock = true;
      work.share = decompose(cloud, rank, nRanks, halo);
      gForceBlock = false;
    }
#else
    work.share = decompose(cloud, rank, nRanks, halo);
#endif
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
    cache.tag = workTag;
    cache.epoch = epoch;
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
#ifdef SEAMS_HAS_LINKCELL
    if (usedCut) {
      const std::vector<int> &assign =
          weighted ? gMeasured[workTag].assign : gCount.assign;
      linkcellSends(cloud, work.share, assign, cache.nbrs, cache.sendAt);
    } else
#endif
    {
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
    // Rings never post boundary values, so the plan stays where the first
    // call of this cut left it. Restoring it on every reuse copied the
    // send lists and then exchanged coordinates the rank already holds.
    if (workTag == kWorkRings) {
      // Filled below from the local frame.
    } else {
    // The value exchange reads one plan. Restoring it here keeps a block
    // call from posting the cell cut's graph, which truncated the message.
    const int cachedDegree = static_cast<int>(cache.nbrs.size());
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
    }
  }
  // Rings read coordinates the rank already holds. The ghost skin is the
  // extra atoms in the neighbour list, not a second copy of those coordinates
  // sent across the distributed graph. Boundary values are a later exchange
  // and rings do not post one.
  if (workTag == kWorkRings) {
    // One filled cloud per cut. Alternating cuts would otherwise default-
    // construct the point vector on every call.
    static Cloud kept[2];
    static std::vector<int> placedRings[2];
    static int resident = -1;
    const int slot = useCell ? 1 : 0;
    if (resident != slot) {
      if (resident >= 0) {
        kept[resident] = std::move(work.cloud);
      }
      work.cloud = std::move(kept[slot]);
      resident = slot;
    }
    const bool fresh = placedRings[slot] != work.share.local ||
                       work.cloud.pts.size() != work.share.local.size();
    if (fresh) {
      placedRings[slot] = work.share.local;
      work.cloud.currentFrame = cloud.currentFrame;
      work.cloud.box = cloud.box;
      work.cloud.boxLow = cloud.boxLow;
      work.cloud.pts.assign(work.share.local.size(), {});
      work.cloud.idIndexMap.clear();
    }
    for (std::size_t k = 0; k < work.share.local.size(); k++) {
      const auto &src = cloud.pts[static_cast<std::size_t>(work.share.local[k])];
      auto &dst = work.cloud.pts[k];
      dst.x = src.x;
      dst.y = src.y;
      dst.z = src.z;
      if (fresh) {
        dst.atomID = src.atomID;
        dst.molID = src.molID;
        dst.type = src.type;
      }
    }
    work.cloud.nop = static_cast<int>(work.cloud.pts.size());
    work.commMs = 0.0;
    return work;
  }
  const std::vector<int> &nbrs = cache.nbrs;
  const int degree = static_cast<int>(nbrs.size());
  MPI_Comm graph = cache.graph;
  // Every rank already holds this frame, so the halo exchange carries the
  // frame index. Coordinates are copied from the local cloud below.
  std::vector<std::vector<int>> bags(static_cast<std::size_t>(degree));
  for (int sidx = 0; sidx < degree; sidx++) {
    auto &bag = bags[static_cast<std::size_t>(sidx)];
    bag.reserve(cache.sendAt[static_cast<std::size_t>(sidx)].size());
    for (const int k : cache.sendAt[static_cast<std::size_t>(sidx)]) {
      bag.push_back(work.share.local[static_cast<std::size_t>(k)]);
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
  std::vector<int> send(static_cast<std::size_t>(std::max(sendN, 1)));
  for (int sidx = 0; sidx < degree; sidx++) {
    std::copy(bags[static_cast<std::size_t>(sidx)].begin(),
              bags[static_cast<std::size_t>(sidx)].end(),
              send.begin() + sendDispl[static_cast<std::size_t>(sidx)]);
  }
  const int atomBytes = static_cast<int>(sizeof(int));
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
  if (const char *split = std::getenv("SEAMS_SPLIT")) {
    if (split[0] == '1' && rank == 0) {
      std::fprintf(stderr, "send atoms %d degree %d\n", sendN, degree);
    }
  }
  MPI_Barrier(comm);
  const auto t0 = std::chrono::steady_clock::now();
  MPI_Neighbor_alltoall(sendBytesPtr, 1, MPI_INT, recvBytes.data(), 1, MPI_INT, graph);
  int recvN = 0;
  for (int sidx = 0; sidx < degree; sidx++) {
    recvByteDispl[static_cast<std::size_t>(sidx)] = recvN;
    recvN += recvBytes[static_cast<std::size_t>(sidx)];
  }
  std::vector<char> recv(static_cast<std::size_t>(std::max(recvN, 1)));
  // One cloud per cut. Alternating block and cell would otherwise
  // default-construct every point, including its bond vector, on each call.
  static Cloud keptCloud[2];
  static std::vector<int> placedCut[2];
  static int residentCut = -1;
  const int cutSlot = useCell ? 1 : 0;
  if (residentCut != cutSlot) {
    if (residentCut >= 0) {
      keptCloud[residentCut] = std::move(work.cloud);
    }
    work.cloud = std::move(keptCloud[cutSlot]);
    residentCut = cutSlot;
  }
  work.cloud.currentFrame = cloud.currentFrame;
  work.cloud.box = cloud.box;
  work.cloud.boxLow = cloud.boxLow;
  if (work.cloud.pts.size() != work.share.local.size()) {
    work.cloud.pts.assign(work.share.local.size(), {});
    work.cloud.idIndexMap.clear();
  }
  const auto copyT0 = std::chrono::steady_clock::now();
  const int nLocalCopy = static_cast<int>(work.share.local.size());
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (omp_get_max_threads() > 1)
#endif
  for (int k = 0; k < nLocalCopy; k++) {
    const int frame = work.share.local[k];
    if (frame < 0 || frame >= cloud.nop) {
      continue;
    }
    const auto &src = cloud.pts[static_cast<std::size_t>(frame)];
    auto &dst = work.cloud.pts[k];
    dst.x = src.x;
    dst.y = src.y;
    dst.z = src.z;
    dst.atomID = src.atomID;
    dst.molID = src.molID;
    dst.type = src.type;
  }
  if (const char *split = std::getenv("SEAMS_SPLIT")) {
    if (split[0] == '1' && rank == 0) {
      std::fprintf(stderr, "copy ms %.3f\n",
                   std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - copyT0)
                       .count());
    }
  }
  work.cloud.nop = static_cast<int>(work.cloud.pts.size());
  MPI_Request payload = MPI_REQUEST_NULL;
  MPI_Ineighbor_alltoallv(send.data(), sendBytesPtr, sendByteDisplPtr, MPI_BYTE, recv.data(),
                          recvBytes.data(), recvByteDispl.data(), MPI_BYTE, graph, &payload);
  // Owned coordinates are already on this rank. Ghost coordinates are the
  // bytes in flight, so a hook can run the owned kernel against that wait.
  if (gHaloHook != nullptr) {
    gHaloHook(gHaloHookCtx);
  }
  MPI_Wait(&payload, MPI_STATUS_IGNORE);
  MPI_Barrier(comm);
  work.commMs = std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - t0)
                    .count();
  (void)graph;

  // Place ghosts by frame index. The local list is sorted, so each
  // received atom is a binary search, and the point is rewritten in place.
  // Two cuts can hold the same number of atoms and still own different ones.
  const bool fresh = placedCut[cutSlot] != work.share.local ||
                     work.cloud.pts.size() != work.share.local.size();
  if (fresh) {
    placedCut[cutSlot] = work.share.local;
  }
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
    int frame = 0;
    std::memcpy(&frame, recv.data() + at, sizeof(int));
    const auto it = std::lower_bound(work.share.local.begin(), work.share.local.end(), frame);
    if (it == work.share.local.end() || *it != frame) {
      misses++;
      continue;
    }
    const std::size_t k = static_cast<std::size_t>(it - work.share.local.begin());
    if (fresh) {
      const auto &src = cloud.pts[static_cast<std::size_t>(frame)];
      auto &dst = work.cloud.pts[k];
      dst.x = src.x;
      dst.y = src.y;
      dst.z = src.z;
      dst.atomID = src.atomID;
      dst.molID = src.molID;
      dst.type = src.type;
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

void captureLocalWork(int workTag, const double *perLocal, int nLocal, MPI_Comm comm) {
#ifdef SEAMS_HAS_LINKCELL
  if (workTag <= 0 || workTag >= 5 || perLocal == nullptr || comm == MPI_COMM_NULL) {
    return;
  }
  MeasuredCut &slot = gMeasured[workTag];
  if (slot.ready || gCount.assign.empty() || gCount.pointCell.empty()) {
    return;
  }
  const Share &share = gatheredWork().share;
  const int ncell = static_cast<int>(gCount.assign.size());
  const int n = std::min(nLocal, static_cast<int>(share.local.size()));
  std::vector<std::uint64_t> mine(static_cast<std::size_t>(ncell), 0);
  for (int k = 0; k < n; k++) {
    if (!share.owned[static_cast<std::size_t>(k)] || !(perLocal[k] > 0.0)) {
      continue;
    }
    const int frame = share.local[static_cast<std::size_t>(k)];
    if (frame < 0 || static_cast<std::size_t>(frame) >= gCount.pointCell.size()) {
      continue;
    }
    const int cell = gCount.pointCell[static_cast<std::size_t>(frame)];
    if (cell < 0 || cell >= ncell) {
      continue;
    }
    const auto quantum = static_cast<std::uint64_t>(perLocal[k] + 0.5);
    mine[static_cast<std::size_t>(cell)] += quantum > 0 ? quantum : 1;
  }
  std::vector<std::uint64_t> sum(mine.size(), 0);
  MPI_Allreduce(mine.data(), sum.data(), ncell, MPI_UINT64_T, MPI_SUM, comm);
  std::uint64_t total = 0;
  for (const auto w : sum) {
    total += w;
  }
  if (total == 0 || gCount.parts < 1) {
    return;
  }
  std::vector<int> assign(static_cast<std::size_t>(ncell), 0);
  if (lc_partition_weights(gCount.nxyz, sum.data(), static_cast<std::size_t>(ncell),
                           static_cast<std::size_t>(gCount.parts), 2, 0.0,
                           assign.data()) != 0) {
    return;
  }
  slot.weight = std::move(sum);
  slot.assign = std::move(assign);
  slot.nop = gCount.nop;
  slot.parts = gCount.parts;
  slot.halo = gCount.halo;
  slot.ready = true;
#else
  (void)workTag;
  (void)perLocal;
  (void)nLocal;
  (void)comm;
#endif
}

const std::vector<std::vector<int>> &distributedRings(const Cloud &cloud, double cutoff,
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
    static thread_local std::vector<std::vector<int>> skewed;
    skewed = rings(cloud, cutoff, maxDepth, rank, nRanks);
    return skewed;
  }
  // ringHops is ringHaloBonds, floor(maxDepth/2). The coordinate skin is
  // that many cutoffs: every primitive ring through an owned atom, and
  // every shortcut of it, stays inside the skin. Owned sources are
  // searched on that halo and mapped back to the frame. The bond ball
  // itself ends near 7.8 A on this crystal, but a skin of that radius
  // made the wall longer, so the gather stays at ringHops cutoffs.
  const Work &work = gatherWork(cloud, ringHops(maxDepth) * cutoff, comm, kWorkRings);
  const bool measure = !gMeasured[kWorkRings].ready;
  std::vector<double> sourceUs;
  if (measure) {
    sourceUs.assign(static_cast<std::size_t>(std::max(work.cloud.nop, 0)), 0.0);
    primitive::setSourceWork(sourceUs.empty() ? nullptr : sourceUs.data());
  }
  const auto &local = ownedRings(work.cloud, cutoff, maxDepth, work.share.owned);
  if (measure) {
    primitive::setSourceWork(nullptr);
    captureLocalWork(kWorkRings, sourceUs.data(), static_cast<int>(sourceUs.size()), comm);
  }
  const auto reduce0 = std::chrono::steady_clock::now();
  thread_local std::vector<std::vector<int>> mapped;
  if (mapped.size() < local.size()) {
    mapped.resize(local.size());
  }
  for (std::size_t r = 0; r < local.size(); r++) {
    auto &dst = mapped[r];
    const auto &src = local[r];
    if (dst.size() != src.size()) {
      dst.resize(src.size());
    }
    for (std::size_t i = 0; i < src.size(); i++) {
      dst[i] = work.share.local[static_cast<std::size_t>(src[i])];
    }
  }
  if (mapped.size() != local.size()) {
    mapped.resize(local.size());
  }
  gRingPhases.reduceMs +=
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - reduce0).count();
  return mapped;
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
