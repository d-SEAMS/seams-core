#include <phase.hpp>
#include <cage_enum.hpp>
#include <franzblau.hpp>
#include <generic.hpp>
#include <neighbours.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

namespace {

constexpr double kXXIa = 20.197;
constexpr double kXXIc = 7.891;
constexpr int kXXIn = 152;
constexpr double kXXIrho = 1.413; // g/cm^3, water mass 18.015

} // namespace

namespace {

bool indexBonded(const std::vector<std::vector<int>> &nList, int i, int j) {
  if (i < 0 || j < 0 || static_cast<std::size_t>(i) >= nList.size()) {
    return false;
  }
  const auto &row = nList[static_cast<std::size_t>(i)];
  return std::find(row.begin(), row.end(), j) != row.end();
}

bool prismStacked(const std::vector<int> &a, const std::vector<int> &b,
                  const std::vector<std::vector<int>> &nList) {
  if (a.size() != 6 || b.size() != 6) {
    return false;
  }
  std::unordered_set<int> vb(b.begin(), b.end());
  for (int v : a) {
    if (vb.count(v) != 0) {
      return false;
    }
  }
  for (int va : a) {
    int hits = 0;
    for (int vb0 : b) {
      if (indexBonded(nList, va, vb0)) {
        ++hits;
      }
    }
    if (hits != 1) {
      return false;
    }
  }
  return true;
}

int channelComponents(const std::vector<std::vector<int>> &six,
                      const std::vector<std::vector<int>> &nList) {
  const int n = static_cast<int>(six.size());
  if (n == 0) {
    return 0;
  }
  std::vector<int> parent(static_cast<std::size_t>(n));
  std::iota(parent.begin(), parent.end(), 0);
  auto find = [&](int x) {
    while (parent[static_cast<std::size_t>(x)] != x) {
      parent[static_cast<std::size_t>(x)] =
          parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];
      x = parent[static_cast<std::size_t>(x)];
    }
    return x;
  };
  auto unite = [&](int x, int y) {
    x = find(x);
    y = find(y);
    if (x != y) {
      parent[static_cast<std::size_t>(y)] = x;
    }
  };
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (prismStacked(six[static_cast<std::size_t>(i)],
                       six[static_cast<std::size_t>(j)], nList)) {
        unite(i, j);
      }
    }
  }
  std::vector<int> size(static_cast<std::size_t>(n), 0);
  for (int i = 0; i < n; ++i) {
    ++size[static_cast<std::size_t>(find(i))];
  }
  int channels = 0;
  for (int s : size) {
    if (s >= 2) {
      ++channels;
    }
  }
  return channels;
}

} // namespace

int phase::openChannelCount(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList) {
  (void)yCloud;
  const auto sig = cage::Signature::parse("512");
  const auto rings7 = primitive::ringNetwork(nList, 7);
  const auto closed = cage::findBySignature(rings7, nList, sig);
  if (!closed.empty()) {
    return 0;
  }
  const auto rings6 = primitive::ringNetwork(nList, 6);
  std::vector<std::vector<int>> six;
  six.reserve(rings6.size());
  for (const auto &r : rings6) {
    if (r.size() == 6) {
      six.push_back(r);
    }
  }
  return channelComponents(six, nList);
}

std::uint64_t phase::protonKey(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    int oxygenType, int hydrogenType) {
  std::unordered_map<int, std::vector<int>> hByMol;
  std::unordered_map<int, int> oByMol;
  for (int i = 0; i < yCloud.nop; i++) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(i)];
    if (p.type == oxygenType) {
      oByMol[p.molID] = i;
    }
    if (p.type == hydrogenType) {
      hByMol[p.molID].push_back(i);
    }
  }
  std::uint64_t key = 1469598103934665603ULL;
  std::vector<int> mols;
  mols.reserve(oByMol.size());
  for (const auto &kv : oByMol) {
    mols.push_back(kv.first);
  }
  std::sort(mols.begin(), mols.end());
  for (int mol : mols) {
    const int oi = oByMol[mol];
    auto hs = hByMol[mol];
    std::sort(hs.begin(), hs.end());
    for (int hi : hs) {
      const auto dr = gen::relDist(yCloud, hi, oi);
      const int bx = static_cast<int>(std::lround(dr[0] * 4.0));
      const int by = static_cast<int>(std::lround(dr[1] * 4.0));
      const int bz = static_cast<int>(std::lround(dr[2] * 4.0));
      const std::uint64_t word =
          (static_cast<std::uint64_t>(static_cast<uint32_t>(bx)) << 32) ^
          (static_cast<std::uint64_t>(static_cast<uint32_t>(by)) << 16) ^
          static_cast<std::uint64_t>(static_cast<uint32_t>(bz));
      key ^= word;
      key *= 1099511628211ULL;
    }
  }
  return key;
}

double phase::hydrogenMSD(
    const molSys::PointCloud<molSys::Point<double>, double> &frame0,
    const molSys::PointCloud<molSys::Point<double>, double> &frame1,
    int hydrogenType) {
  std::unordered_map<int, int> idx0;
  for (int i = 0; i < frame0.nop; ++i) {
    const auto &p = frame0.pts[static_cast<std::size_t>(i)];
    if (p.type == hydrogenType) {
      idx0[p.atomID] = i;
    }
  }
  double acc = 0.0;
  int n = 0;
  for (const auto &p : frame1.pts) {
    if (p.type != hydrogenType) {
      continue;
    }
    const auto it = idx0.find(p.atomID);
    if (it == idx0.end()) {
      continue;
    }
    const auto dr = gen::relDistFromPoint(frame0, it->second, p.x, p.y, p.z);
    acc += dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
    ++n;
  }
  return n > 0 ? acc / static_cast<double>(n)
               : std::numeric_limits<double>::quiet_NaN();
}

phase::IceXXIHit phase::iceXXILibrary(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  IceXXIHit hit;
  hit.nSites = yCloud.nop;
  if (yCloud.box.size() >= 3) {
    hit.a = 0.5 * (yCloud.box[0] + yCloud.box[1]);
    hit.c = yCloud.box[2];
  }
  const double vol = (yCloud.box.size() >= 3)
                         ? yCloud.box[0] * yCloud.box[1] * yCloud.box[2]
                         : 0.0;
  const double mass = static_cast<double>(yCloud.nop) * 18.015 / 6.02214076e23;
  hit.density = vol > 0.0 ? (mass / (vol * 1e-24)) : 0.0;
  const bool nOk = yCloud.nop == kXXIn;
  const bool aOk = std::fabs(hit.a - kXXIa) < 0.4;
  const bool cOk = std::fabs(hit.c - kXXIc) < 0.3;
  const bool rhoOk = std::fabs(hit.density - kXXIrho) < 0.08;
  // Ice I uses 3.5 A; ice XXI at 1.413 g/cm3 has a 3.0 A first shell.
  // Four-nearest neighbours are the tetrahedral graph in both cells.
  auto nList = nneigh::kNearestNeighbourList(yCloud, 4, 3.5, 1, true);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  double coord = 0.0;
  for (const auto &row : nList) {
    if (!row.empty()) {
      coord += static_cast<double>(row.size() - 1);
    }
  }
  hit.meanCoord = nList.empty() ? 0.0 : coord / static_cast<double>(nList.size());
  const auto rings = primitive::ringNetwork(nList, 6);
  for (const auto &r : rings) {
    hit.nSix += static_cast<int>(r.size() == 6);
  }
  const bool tetra = hit.meanCoord >= 3.5 && hit.meanCoord <= 4.5;
  // The Lee cell has ~112 primitive six-rings on the 4-NN graph.
  // A cubic packing of the same box can form one accidental six-ring.
  const bool ringsOk = hit.nSix >= 50;
  hit.match = nOk && aOk && cOk && rhoOk && tetra && ringsOk;
  return hit;
}

std::vector<double> phase::localDensity(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    double rcut) {
  std::vector<double> rho(static_cast<std::size_t>(yCloud.nop), 0.0);
  const double r2 = rcut * rcut;
  const double vol = (4.0 / 3.0) * 3.14159265358979323846 * rcut * rcut * rcut;
  for (int i = 0; i < yCloud.nop; i++) {
    int n = 0;
    for (int j = 0; j < yCloud.nop; j++) {
      if (i == j) {
        continue;
      }
      const auto dr = gen::relDist(yCloud, j, i);
      if (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] < r2) {
        ++n;
      }
    }
    rho[static_cast<std::size_t>(i)] = static_cast<double>(n) / vol;
  }
  return rho;
}

double phase::frameDensity(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  const double vol = (yCloud.box.size() >= 3)
                         ? yCloud.box[0] * yCloud.box[1] * yCloud.box[2]
                         : 0.0;
  const double mass = static_cast<double>(yCloud.nop) * 18.015 / 6.02214076e23;
  return vol > 0.0 ? (mass / (vol * 1e-24)) : 0.0;
}

phase::GlassKind phase::glassFromDensity(double rho, double iceMax,
                                         double ldaMax, double mdaMin) {
  if (rho < iceMax) {
    return GlassKind::ice;
  }
  if (rho < ldaMax) {
    return GlassKind::lda;
  }
  if (rho < mdaMin) {
    return GlassKind::mda;
  }
  return GlassKind::hda;
}

phase::GlassKind phase::glassFromDensity(double rho, double iceMax,
                                         double mdaMin) {
  return glassFromDensity(rho, iceMax, 0.5 * (iceMax + mdaMin), mdaMin);
}

phase::GlassKind phase::glassFromDensity(double rho) {
  return glassFromDensity(rho, 0.028, 0.040, 0.070);
}
