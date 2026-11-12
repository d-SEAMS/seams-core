//-----------------------------------------------------------------------------------
// d-SEAMS
// Copyright (c) d-SEAMS developers
//-----------------------------------------------------------------------------------

#include <diamond_cna.hpp>
#include <generic.hpp>

#include <algorithm>
#include <cmath>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

int diamondCna(const double *sx, const double *sy, const double *sz) {
  double rsum = 0.0;
  for (int i = 0; i < 12; ++i) {
    rsum += std::sqrt(sx[i] * sx[i] + sy[i] * sy[i] + sz[i] * sz[i]);
  }
  const double cut = (rsum / 12.0) * 1.30;
  const double cut2 = cut * cut;
  int n421 = 0;
  int n422 = 0;
  int nBad = 0;
  for (int j = 0; j < 12; ++j) {
    int common[12];
    int nc = 0;
    for (int k = 0; k < 12; ++k) {
      if (k == j) {
        continue;
      }
      const double dx = sx[j] - sx[k];
      const double dy = sy[j] - sy[k];
      const double dz = sz[j] - sz[k];
      if (dx * dx + dy * dy + dz * dz < cut2) {
        if (nc < 12) {
          common[nc] = k;
        }
        ++nc;
      }
    }
    if (nc != 4) {
      ++nBad;
      continue;
    }
    int nBonds = 0;
    int deg[4] = {};
    for (int a = 0; a < 4; ++a) {
      for (int b = a + 1; b < 4; ++b) {
        const int ia = common[a];
        const int ib = common[b];
        const double dx = sx[ia] - sx[ib];
        const double dy = sy[ia] - sy[ib];
        const double dz = sz[ia] - sz[ib];
        if (dx * dx + dy * dy + dz * dz < cut2) {
          ++nBonds;
          ++deg[a];
          ++deg[b];
        }
      }
    }
    if (nBonds != 2) {
      ++nBad;
      continue;
    }
    int chain = 1;
    for (int a = 0; a < 4; ++a) {
      if (deg[a] >= 2) {
        chain = 2;
      }
    }
    if (chain == 1) {
      ++n421;
    } else {
      ++n422;
    }
  }
  if (n421 == 12) {
    return 1;
  }
  if (n421 == 6 && n422 == 6) {
    return 0;
  }
  if (n421 + n422 >= 4 && nBad <= 6) {
    return -2;
  }
  return -1;
}

} // namespace

namespace symm {

void diamondFromCsr(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                    const std::vector<int> &offsets, const std::vector<int> &cols,
                    std::vector<int> &kind) {
  const int n = yCloud.nop;
  kind.assign(static_cast<std::size_t>(std::max(n, 0)), -1);
  if (n <= 0 || static_cast<int>(offsets.size()) < n + 1) {
    return;
  }
  const gen::FracBox frame = gen::makeFracBox(yCloud);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= 256)
#endif
  for (int i = 0; i < n; ++i) {
    struct Bond {
      int j;
      double d2, x, y, z;
    };
    Bond row[16];
    int nb = 0;
    const auto &pi = yCloud.pts[static_cast<std::size_t>(i)];
    const int a = offsets[static_cast<std::size_t>(i)];
    const int b = offsets[static_cast<std::size_t>(i) + 1];
    for (int k = a; k < b && nb < 16; ++k) {
      const int j = cols[static_cast<std::size_t>(k)];
      if (j < 0 || j >= n || j == i) {
        continue;
      }
      const auto &pj = yCloud.pts[static_cast<std::size_t>(j)];
      const auto dr = gen::fracDelta(frame, pi.x, pi.y, pi.z, pj.x, pj.y, pj.z);
      const double d2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
      if (!(d2 > 2.15 * 2.15)) {
        continue;
      }
      row[nb++] = {j, d2, dr[0], dr[1], dr[2]};
    }
    if (nb < 4) {
      continue;
    }
    std::sort(row, row + nb, [](const Bond &p, const Bond &q) { return p.d2 < q.d2; });
    double u[4][3];
    double acc = 0.0;
    bool tetra = true;
    for (int k = 0; k < 4; ++k) {
      const double nrm = std::sqrt(row[k].d2);
      if (!(nrm > 1e-8)) {
        tetra = false;
        break;
      }
      u[k][0] = row[k].x / nrm;
      u[k][1] = row[k].y / nrm;
      u[k][2] = row[k].z / nrm;
    }
    if (!tetra) {
      continue;
    }
    for (int p = 0; p < 4; ++p) {
      for (int q = p + 1; q < 4; ++q) {
        const double c = u[p][0] * u[q][0] + u[p][1] * u[q][1] + u[p][2] * u[q][2];
        const double d = c + (1.0 / 3.0);
        acc += d * d;
      }
    }
    if ((1.0 - 0.375 * acc) < 0.88) {
      continue;
    }
    double fx[12], fy[12], fz[12];
    int nf = 0;
    bool closed = true;
    for (int k = 0; k < 4 && closed; ++k) {
      const int f = row[k].j;
      const int fa = offsets[static_cast<std::size_t>(f)];
      const int fb = offsets[static_cast<std::size_t>(f) + 1];
      int taken = 0;
      for (int t = fa; t < fb && taken < 3; ++t) {
        const int j = cols[static_cast<std::size_t>(t)];
        if (j == i || j == row[0].j || j == row[1].j || j == row[2].j || j == row[3].j) {
          continue;
        }
        const auto &pj = yCloud.pts[static_cast<std::size_t>(j)];
        const auto &pf = yCloud.pts[static_cast<std::size_t>(f)];
        const auto arm = gen::fracDelta(frame, pf.x, pf.y, pf.z, pj.x, pj.y, pj.z);
        const double x = row[k].x + arm[0];
        const double y = row[k].y + arm[1];
        const double z = row[k].z + arm[2];
        const double origin = x * x + y * y + z * z;
        if (origin < 0.25) {
          continue;
        }
        bool isFirst = false;
        for (int q = 0; q < 4; ++q) {
          const double dx = x - row[q].x;
          const double dy = y - row[q].y;
          const double dz = z - row[q].z;
          if (dx * dx + dy * dy + dz * dz < 0.25) {
            isFirst = true;
          }
        }
        if (isFirst || nf >= 12) {
          if (nf >= 12) {
            closed = false;
          }
          continue;
        }
        fx[nf] = x;
        fy[nf] = y;
        fz[nf] = z;
        ++nf;
        ++taken;
      }
      if (taken < 3) {
        closed = false;
      }
    }
    if (!(closed && nf == 12)) {
      continue;
    }
    const int cna = diamondCna(fx, fy, fz);
    if (cna >= 0) {
      kind[static_cast<std::size_t>(i)] = cna;
    }
  }
}

} // namespace symm
