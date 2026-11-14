// Five-seed Ic/Ih/null sweep for Rodger F4 and Zeron q3/q12.
// Same lattices as scripts/sota_compare.py (512 cubic, 360 hexagonal,
// 360 null). Decision rules are fit on the ideal lattices.
#include <bop.hpp>
#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <order_parameter.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double kBond = 2.75;
constexpr std::uint64_t kSeed = 88172645463325252ULL;
const double kSigmas[] = {0.20, 0.30, 0.35, 0.40};

double wrapLen(double x, double L) {
  x = std::fmod(x, L);
  if (x < 0.0) {
    x += L;
  }
  return x;
}

molSys::PointCloud<molSys::Point<double>, double>
cloudFromPos(const std::vector<std::array<double, 3>> &pos,
             const std::array<double, 3> &box) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {box[0], box[1], box[2]};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  for (std::size_t i = 0; i < pos.size(); i++) {
    molSys::Point<double> p;
    p.type = 1;
    p.atomID = static_cast<int>(i) + 1;
    p.molID = p.atomID;
    p.x = pos[i][0];
    p.y = pos[i][1];
    p.z = pos[i][2];
    cloud.pts.push_back(p);
    cloud.idIndexMap[p.atomID] = static_cast<int>(i);
  }
  cloud.nop = static_cast<int>(pos.size());
  return cloud;
}

std::vector<std::array<double, 3>>
cubicDiamond(int reps, std::array<double, 3> &box) {
  const double a = 4.0 * kBond / std::sqrt(3.0);
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  std::vector<std::array<double, 3>> pos;
  for (int i = 0; i < reps; i++) {
    for (int j = 0; j < reps; j++) {
      for (int k = 0; k < reps; k++) {
        for (int b = 0; b < 4; b++) {
          for (double s : {0.0, 0.25}) {
            pos.push_back({wrapLen((i + fcc[b][0] + s) * a, reps * a),
                           wrapLen((j + fcc[b][1] + s) * a, reps * a),
                           wrapLen((k + fcc[b][2] + s) * a, reps * a)});
          }
        }
      }
    }
  }
  box = {reps * a, reps * a, reps * a};
  return pos;
}

std::vector<std::array<double, 3>>
lonsdaleite(int nx, int ny, int nz, std::array<double, 3> &box) {
  const double a = kBond / ((3.0 / 8.0) * std::sqrt(8.0 / 3.0));
  const double c = a * std::sqrt(8.0 / 3.0);
  const double hex[4][3] = {{1.0 / 3.0, 2.0 / 3.0, 0.0},
                            {2.0 / 3.0, 1.0 / 3.0, 0.5},
                            {1.0 / 3.0, 2.0 / 3.0, 3.0 / 8.0},
                            {2.0 / 3.0, 1.0 / 3.0, 7.0 / 8.0}};
  const double X = a;
  const double Y = a * std::sqrt(3.0);
  std::vector<std::array<double, 3>> ortho;
  for (const auto &h : hex) {
    const double x = (h[0] - 0.5 * h[1]) * a;
    const double y = h[1] * (std::sqrt(3.0) / 2.0) * a;
    for (auto [sx, sy] : {std::pair{0.0, 0.0}, std::pair{X / 2.0, Y / 2.0}}) {
      ortho.push_back({wrapLen(x + sx, X), wrapLen(y + sy, Y), h[2] * c});
    }
  }
  box = {nx * X, ny * Y, nz * c};
  std::vector<std::array<double, 3>> pos;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        for (const auto &o : ortho) {
          pos.push_back({wrapLen(i * X + o[0], box[0]),
                         wrapLen(j * Y + o[1], box[1]),
                         wrapLen(k * c + o[2], box[2])});
        }
      }
    }
  }
  return pos;
}

std::vector<std::array<double, 3>>
denseNull(int n, const std::array<double, 3> &box, std::mt19937_64 &rng) {
  std::uniform_real_distribution<double> ux(0.0, box[0]);
  std::uniform_real_distribution<double> uy(0.0, box[1]);
  std::uniform_real_distribution<double> uz(0.0, box[2]);
  std::vector<std::array<double, 3>> pos;
  const double cell2 = 2.3 * 2.3;
  while (static_cast<int>(pos.size()) < n) {
    const std::array<double, 3> cand{ux(rng), uy(rng), uz(rng)};
    bool ok = true;
    const int start = std::max(0, static_cast<int>(pos.size()) - 400);
    for (int i = start; i < static_cast<int>(pos.size()); i++) {
      double d2 = 0.0;
      for (std::size_t k = 0; k < 3; k++) {
        double d = std::abs(cand[k] - pos[static_cast<std::size_t>(i)][k]);
        d = std::min(d, box[k] - d);
        d2 += d * d;
      }
      if (d2 < cell2) {
        ok = false;
        break;
      }
    }
    if (ok) {
      pos.push_back(cand);
    }
  }
  return pos;
}

std::vector<std::array<double, 3>>
jitter(const std::vector<std::array<double, 3>> &pos,
       const std::array<double, 3> &box, double sigma, std::uint64_t seed) {
  if (sigma <= 0.0) {
    return pos;
  }
  std::mt19937_64 rng(seed);
  std::normal_distribution<double> g(0.0, sigma);
  auto out = pos;
  for (auto &p : out) {
    for (std::size_t k = 0; k < 3; k++) {
      p[k] = wrapLen(p[k] + g(rng), box[k]);
    }
  }
  return out;
}

struct QPair {
  double q3 = 0.0;
  double q12 = 0.0;
};

std::vector<QPair> zeronOnOxygens(const std::vector<std::array<double, 3>> &pos,
                                  const std::array<double, 3> &box) {
  auto cloud = cloudFromPos(pos, box);
  auto nList = nneigh::kNearestNeighbourList(cloud, 4, 5.0, 1, true);
  const auto q3 = chill::steinhardtQl(cloud, nList, 3);
  const auto q12 = chill::steinhardtQl(cloud, nList, 12);
  std::vector<QPair> out(pos.size());
  for (std::size_t i = 0; i < pos.size(); i++) {
    out[i].q3 = q3.qlBar[i];
    out[i].q12 = q12.qlBar[i];
  }
  return out;
}

QPair meanPair(const std::vector<QPair> &v) {
  QPair m;
  int n = 0;
  for (const auto &p : v) {
    if (std::isfinite(p.q3) && std::isfinite(p.q12)) {
      m.q3 += p.q3;
      m.q12 += p.q12;
      ++n;
    }
  }
  if (n > 0) {
    m.q3 /= static_cast<double>(n);
    m.q12 /= static_cast<double>(n);
  }
  return m;
}

double dist2(QPair a, QPair b) {
  const double dq3 = a.q3 - b.q3;
  const double dq12 = a.q12 - b.q12;
  return dq3 * dq3 + dq12 * dq12;
}

void report(const std::string &method, const std::string &system, double sigma,
            int seed, double crystal, const std::string &note) {
  std::cout << "method=" << method << " system=" << system << " sigma="
            << std::fixed << std::setprecision(2) << sigma << " seed=" << seed
            << " acc=nan cubic=nan hex=nan crystal=" << std::setprecision(3)
            << crystal << " note=" << note << "\n";
}

} // namespace

int main() {
  std::array<double, 3> icBox{};
  std::array<double, 3> ihBox{};
  const auto ic0 = cubicDiamond(4, icBox);
  const auto ih0 = lonsdaleite(5, 3, 3, ihBox);
  std::mt19937_64 nullRng(kSeed);
  const auto null0 = denseNull(360, icBox, nullRng);
  std::cout << "# ic n=" << ic0.size() << " ih n=" << ih0.size()
            << " null n=" << null0.size() << "\n";

  std::cout << "# rodger-f4 undefined on these oxygen-only mW lattices\n";

  const auto zIc0 = zeronOnOxygens(ic0, icBox);
  const auto zIh0 = zeronOnOxygens(ih0, ihBox);
  const auto zN0 = zeronOnOxygens(null0, icBox);
  const QPair cIc = meanPair(zIc0);
  const QPair cIh = meanPair(zIh0);
  const QPair cW = meanPair(zN0);
  const double q3Thr = 0.5 * (0.5 * (cIc.q3 + cIh.q3) + cW.q3);
  std::cout << "# zeron_ic=(" << cIc.q3 << "," << cIc.q12 << ") ih=(" << cIh.q3
            << "," << cIh.q12 << ") null=(" << cW.q3 << "," << cW.q12
            << ") q3_thr=" << q3Thr << "\n";

  struct Sys {
    const char *name;
    const std::vector<std::array<double, 3>> *pos;
    const std::array<double, 3> *box;
  };
  const Sys systems[] = {{"ic", &ic0, &icBox}, {"ih", &ih0, &ihBox},
                         {"null", &null0, &icBox}};

  for (const auto &sys : systems) {
    for (int seed = 0; seed < 5; seed++) {
      for (double sigma : kSigmas) {
        const std::uint64_t js =
            kSeed + static_cast<std::uint64_t>(sigma * 1000.0) +
            7919ULL * static_cast<std::uint64_t>(seed);
        const auto pos = jitter(*sys.pos, *sys.box, sigma, js);
        report("rodger-f4", sys.name, sigma, seed, std::numeric_limits<double>::quiet_NaN(),
               "undefined-on-mW");

        const auto z = zeronOnOxygens(pos, *sys.box);
        int nQ3 = 0;
        int nCentIce = 0;
        for (const auto &p : z) {
          if (!std::isfinite(p.q3) || !std::isfinite(p.q12)) {
            continue;
          }
          if (p.q3 > q3Thr) {
            ++nQ3;
          }
          const double dIc = dist2(p, cIc);
          const double dIh = dist2(p, cIh);
          const double dW = dist2(p, cW);
          if (std::min(dIc, dIh) < dW) {
            ++nCentIce;
          }
        }
        const double den = static_cast<double>(z.size());
        report("zeron-q3", sys.name, sigma, seed, static_cast<double>(nQ3) / den,
               "ice-vs-not");
        report("zeron-q3q12", sys.name, sigma, seed,
               static_cast<double>(nCentIce) / den, "centroid-ice-vs-not");
      }
    }
  }
  return 0;
}
