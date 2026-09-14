#include <chill_offload.hpp>

#include <bop.hpp>
#include <generic.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <utility>
#include <vector>

#ifdef SEAMS_HAS_OFFLOAD
#include <omp.h>
#endif

namespace {

constexpr int kL = 3;
constexpr int kM = 7;
constexpr int kNb = 4;

bool envOffload() {
  const char *v = std::getenv("SEAMS_OFFLOAD");
  return v != nullptr && v[0] != '\0' && !(v[0] == '0' && v[1] == '\0');
}

void y3(double dx, double dy, double dz, double *re, double *im) {
  const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
  const double phi = std::atan2(dx, dy);
  const double theta = (r > 0.0) ? std::acos(dz / r) : 0.0;
  const double s = std::sin(theta);
  const double c = std::cos(theta);
  const double s2 = s * s;
  const double s3 = s2 * s;
  const double c2 = c * c;
  const double c3 = c2 * c;
  constexpr double pi = 3.14159265358979323846;
  const double n0 = 0.25 * std::sqrt(7.0 / pi);
  const double n1 = 0.125 * std::sqrt(21.0 / pi);
  const double n2 = 0.25 * std::sqrt(105.0 / (2.0 * pi));
  const double n3 = 0.125 * std::sqrt(35.0 / pi);
  const double amp0 = n0 * (5.0 * c3 - 3.0 * c);
  const double amp1 = n1 * s * (5.0 * c2 - 1.0);
  const double amp2 = n2 * s2 * c;
  const double amp3 = n3 * s3;
  const double cp = std::cos(phi);
  const double sp = std::sin(phi);
  const double c2p = cp * cp - sp * sp;
  const double s2p = 2.0 * cp * sp;
  const double c3p = cp * c2p - sp * s2p;
  const double s3p = sp * c2p + cp * s2p;
  // m = -3 .. +3 at slots 0..6. harmonicPair: Y_{-m} = amp * conj(phase^m),
  // Y_{+m} = ((m even) ? amp : -amp) * phase^m
  re[3] = amp0;
  im[3] = 0.0;
  re[2] = amp1 * cp;
  im[2] = -amp1 * sp;
  re[4] = -amp1 * cp;
  im[4] = -amp1 * sp;
  re[1] = amp2 * c2p;
  im[1] = -amp2 * s2p;
  re[5] = amp2 * c2p;
  im[5] = amp2 * s2p;
  re[0] = amp3 * c3p;
  im[0] = -amp3 * s3p;
  re[6] = -amp3 * c3p;
  im[6] = -amp3 * s3p;
}

int classifyAtom(int nstag, int neclip, int nb, const int *neigh,
                 const int *nstagN) {
  if (nb != 4) {
    return static_cast<int>(molSys::atom_state_type::water);
  }
  if (neclip == 0 && nstag == 4) {
    return static_cast<int>(molSys::atom_state_type::cubic);
  }
  if (neclip == 1 && nstag == 3) {
    return static_cast<int>(molSys::atom_state_type::hexagonal);
  }
  bool inter = false;
  if (nstag == 2) {
    for (int t = 0; t < 4; t++) {
      const int j = neigh[t];
      if (j >= 0 && nstagN[j] > 2) {
        inter = true;
      }
    }
  }
  if (nstag == 3 && neclip == 0) {
    for (int t = 0; t < 4; t++) {
      const int j = neigh[t];
      if (j >= 0 && nstagN[j] > 1) {
        inter = true;
      }
    }
  }
  if (inter) {
    return static_cast<int>(molSys::atom_state_type::interfacial);
  }
  if (neclip == 4 && nstag == 0) {
    return static_cast<int>(molSys::atom_state_type::clathrate);
  }
  if (neclip == 3) {
    return static_cast<int>(molSys::atom_state_type::interClathrate);
  }
  return static_cast<int>(molSys::atom_state_type::water);
}

} // namespace

bool chill::preferOffload() { return envOffload(); }

chill::ChillPlusResult chill::hostChillPlus(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList) {
  chill::getCorrelPlus(yCloud, nList, false);
  chill::getIceTypePlusNoPrint(yCloud, nList, false);
  ChillPlusResult out;
  out.usedDevice = false;
  out.iceType.resize(static_cast<std::size_t>(yCloud.nop));
  for (int i = 0; i < yCloud.nop; i++) {
    out.iceType[static_cast<std::size_t>(i)] =
        static_cast<int>(yCloud.pts[static_cast<std::size_t>(i)].iceType);
  }
  return out;
}

chill::ChillPlusResult chill::specializedChillPlus(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList) {
  ChillPlusResult out;
  const int n = yCloud.nop;
  out.iceType.assign(static_cast<std::size_t>(n),
                     static_cast<int>(molSys::atom_state_type::unclassified));
  if (n <= 0) {
    return out;
  }
  std::vector<int> neigh(static_cast<std::size_t>(n * kNb), -1);
  std::vector<int> nb(static_cast<std::size_t>(n), 0);
  std::vector<double> xyz(static_cast<std::size_t>(n * 3));
  for (int i = 0; i < n; i++) {
    xyz[static_cast<std::size_t>(3 * i)] = yCloud.pts[static_cast<std::size_t>(i)].x;
    xyz[static_cast<std::size_t>(3 * i + 1)] = yCloud.pts[static_cast<std::size_t>(i)].y;
    xyz[static_cast<std::size_t>(3 * i + 2)] = yCloud.pts[static_cast<std::size_t>(i)].z;
    std::vector<std::pair<double, int>> cand;
    if (static_cast<std::size_t>(i) < nList.size()) {
      const auto &row = nList[static_cast<std::size_t>(i)];
      for (std::size_t t = 1; t < row.size(); t++) {
        const auto it = yCloud.idIndexMap.find(row[t]);
        if (it == yCloud.idIndexMap.end()) {
          continue;
        }
        cand.emplace_back(gen::periodicDistSq(yCloud, i, it->second),
                          it->second);
      }
    }
    const std::size_t keep = std::min(static_cast<std::size_t>(kNb), cand.size());
    std::partial_sort(cand.begin(), cand.begin() + keep, cand.end());
    for (std::size_t t = 0; t < keep; t++) {
      neigh[static_cast<std::size_t>(i * kNb + t)] = cand[t].second;
    }
    nb[static_cast<std::size_t>(i)] = static_cast<int>(keep);
  }
  const double box[3] = {yCloud.box.empty() ? 0.0 : yCloud.box[0],
                         yCloud.box.size() > 1 ? yCloud.box[1] : 0.0,
                         yCloud.box.size() > 2 ? yCloud.box[2] : 0.0};

  std::vector<double> qre(static_cast<std::size_t>(n * kM), 0.0);
  std::vector<double> qim(static_cast<std::size_t>(n * kM), 0.0);
  std::vector<int> nstag(static_cast<std::size_t>(n), 0);
  std::vector<int> neclip(static_cast<std::size_t>(n), 0);

  bool onDevice = false;
#ifdef SEAMS_HAS_OFFLOAD
  if (envOffload() && omp_get_num_devices() > 0) {
    onDevice = true;
  }
#endif
  out.usedDevice = onDevice;

  auto mic = [&](double d, double L) {
    if (L <= 0.0) {
      return d;
    }
    return d - L * std::round(d / L);
  };

#ifdef SEAMS_HAS_OFFLOAD
#pragma omp target teams distribute parallel for if (onDevice) \
    map(to : xyz[0 : n * 3], neigh[0 : n * kNb], nb[0 : n], box[0 : 3]) \
    map(from : qre[0 : n * kM], qim[0 : n * kM])
#endif
#ifndef SEAMS_HAS_OFFLOAD
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for
#endif
#endif
  for (int i = 0; i < n; i++) {
    double accRe[kM] = {};
    double accIm[kM] = {};
    const int nbi = nb[static_cast<std::size_t>(i)];
    for (int t = 0; t < nbi; t++) {
      const int j = neigh[static_cast<std::size_t>(i * kNb + t)];
      if (j < 0 || j >= n) {
        continue;
      }
      const double dx = mic(xyz[static_cast<std::size_t>(3 * j)] -
                                xyz[static_cast<std::size_t>(3 * i)],
                            box[0]);
      const double dy = mic(xyz[static_cast<std::size_t>(3 * j + 1)] -
                                xyz[static_cast<std::size_t>(3 * i + 1)],
                            box[1]);
      const double dz = mic(xyz[static_cast<std::size_t>(3 * j + 2)] -
                                xyz[static_cast<std::size_t>(3 * i + 2)],
                            box[2]);
      double re[kM], im[kM];
      y3(dx, dy, dz, re, im);
      for (int m = 0; m < kM; m++) {
        accRe[m] += re[m];
        accIm[m] += im[m];
      }
    }
    const double inv = (nbi > 0) ? 1.0 / static_cast<double>(nbi) : 0.0;
    for (int m = 0; m < kM; m++) {
      qre[static_cast<std::size_t>(i * kM + m)] = accRe[m] * inv;
      qim[static_cast<std::size_t>(i * kM + m)] = accIm[m] * inv;
    }
  }

#ifdef SEAMS_HAS_OFFLOAD
#pragma omp target teams distribute parallel for if (onDevice) \
    map(to : qre[0 : n * kM], qim[0 : n * kM], neigh[0 : n * kNb], nb[0 : n]) \
    map(from : nstag[0 : n], neclip[0 : n])
#endif
#ifndef SEAMS_HAS_OFFLOAD
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for
#endif
#endif
  for (int i = 0; i < n; i++) {
    int ns = 0;
    int ne = 0;
    const int nbi = nb[static_cast<std::size_t>(i)];
    for (int t = 0; t < nbi; t++) {
      const int j = neigh[static_cast<std::size_t>(i * kNb + t)];
      if (j < 0 || j >= n) {
        continue;
      }
      double dotR = 0.0, iN = 0.0, jN = 0.0;
      for (int m = 0; m < kM; m++) {
        const double qiR = qre[static_cast<std::size_t>(i * kM + m)];
        const double qiI = qim[static_cast<std::size_t>(i * kM + m)];
        const double qjR = qre[static_cast<std::size_t>(j * kM + m)];
        const double qjI = qim[static_cast<std::size_t>(j * kM + m)];
        dotR += qiR * qjR + qiI * qjI;
        iN += qiR * qiR + qiI * qiI;
        jN += qjR * qjR + qjI * qjI;
      }
      const double den = std::sqrt(iN * jN);
      const double cij = (den > 0.0) ? (dotR / den) : 0.0;
      if (cij <= -0.8) {
        ++ns;
      } else if (cij >= -0.35 && cij <= 0.25) {
        ++ne;
      }
    }
    nstag[static_cast<std::size_t>(i)] = ns;
    neclip[static_cast<std::size_t>(i)] = ne;
  }

  for (int i = 0; i < n; i++) {
    out.iceType[static_cast<std::size_t>(i)] = classifyAtom(
        nstag[static_cast<std::size_t>(i)], neclip[static_cast<std::size_t>(i)],
        nb[static_cast<std::size_t>(i)],
        neigh.data() + static_cast<std::size_t>(i * kNb), nstag.data());
    yCloud.pts[static_cast<std::size_t>(i)].iceType =
        static_cast<molSys::atom_state_type>(
            out.iceType[static_cast<std::size_t>(i)]);
  }
  return out;
}

chill::ChillPlusResult chill::chillPlus(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList) {
  if (preferOffload()) {
    return specializedChillPlus(yCloud, nList);
  }
  return hostChillPlus(yCloud, nList);
}
