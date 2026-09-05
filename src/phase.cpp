#include <phase.hpp>
#include <cage_enum.hpp>
#include <franzblau.hpp>
#include <generic.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>

namespace {

constexpr double kXXIa = 20.197;
constexpr double kXXIc = 7.891;
constexpr int kXXIn = 152;
constexpr double kXXIrho = 1.413; // g/cm^3, water mass 18.015

} // namespace

int phase::openChannelCount(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList) {
  const auto sig = cage::Signature::parse("512");
  const auto rings = primitive::ringNetwork(nList, 7);
  const auto closed = cage::findBySignature(rings, nList, sig);
  if (!closed.empty()) {
    return 0;
  }
  int six = 0;
  for (const auto &r : rings) {
    six += static_cast<int>(r.size() == 6);
  }
  return six > 0 ? six : 0;
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
  std::unordered_map<int, std::array<double, 3>> a;
  for (const auto &p : frame0.pts) {
    if (p.type == hydrogenType) {
      a[p.atomID] = {p.x, p.y, p.z};
    }
  }
  double acc = 0.0;
  int n = 0;
  for (const auto &p : frame1.pts) {
    if (p.type != hydrogenType) {
      continue;
    }
    const auto it = a.find(p.atomID);
    if (it == a.end()) {
      continue;
    }
    const double dx = p.x - it->second[0];
    const double dy = p.y - it->second[1];
    const double dz = p.z - it->second[2];
    acc += dx * dx + dy * dy + dz * dz;
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
  hit.match = nOk && aOk && cOk && rhoOk;
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

phase::GlassKind phase::glassFromDensity(double rho, double iceMax,
                                         double mdaMin) {
  if (rho < iceMax) {
    return GlassKind::ice;
  }
  if (rho >= mdaMin) {
    return GlassKind::hda;
  }
  return GlassKind::mda;
}
