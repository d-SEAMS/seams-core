//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
//-----------------------------------------------------------------------------------

#include <neighbours.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

// The linked cell this tree was built against does not always export
// the thread-width hooks. A strong definition in that library wins.
extern "C" __attribute__((weak)) void lc_set_threads(int n) { (void)n; }
extern "C" __attribute__((weak)) int lc_num_threads(void) {
#ifdef SEAMS_HAS_OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

namespace {

constexpr int kMaxCells = 262144;

struct UninitD {
  double v;
  UninitD() noexcept {}
};

int teamRank() {
#ifdef SEAMS_HAS_OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

int teamSize() {
#ifdef SEAMS_HAS_OPENMP
  return omp_get_num_threads();
#else
  return 1;
#endif
}

double minImage2(double dx, double dy, double dz, double lx, double ly,
                 double lz) {
  dx -= lx * std::floor(dx / lx + 0.5);
  dy -= ly * std::floor(dy / ly + 0.5);
  dz -= lz * std::floor(dz / lz + 0.5);
  return dx * dx + dy * dy + dz * dz;
}

void sortRow(int *a, int n) {
  for (int i = 1; i < n; ++i) {
    const int v = a[i];
    int j = i;
    while (j > 0 && a[j - 1] > v) {
      a[j] = a[j - 1];
      --j;
    }
    a[j] = v;
  }
}

struct Scratch {
  std::vector<UninitD> px, py, pz;
  std::vector<UninitD> xs, ys, zs;
  std::vector<int> cellOf;
  std::vector<int> ids;
  std::vector<int> cellStart;
  std::vector<int> hist;
  std::vector<int> counts;
  std::vector<int> cursor;
};

Scratch &scratch() {
  static Scratch s;
  return s;
}

bool profileOn() {
  const char *e = std::getenv("SEAMS_NEIGH_PROFILE");
  return e != nullptr && e[0] != '\0' && e[0] != '0';
}

struct Phase {
  double bin = 0;
  double count = 0;
  double prefix = 0;
  double fill = 0;
  double rows = 0;
  int threads = 1;
  int n = 0;
  int pairs = 0;
};

Phase &phase() {
  static Phase p;
  return p;
}

using Clock = std::chrono::steady_clock;

double msSince(Clock::time_point t0) {
  return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

void publishProfile() {
  if (!profileOn()) {
    return;
  }
  const Phase &p = phase();
  std::fprintf(stderr,
               "neigh_profile n=%d threads=%d pairs=%d bin_ms=%.4f "
               "count_ms=%.4f prefix_ms=%.4f fill_ms=%.4f rows_ms=%.4f\n",
               p.n, p.threads, p.pairs, p.bin, p.count, p.prefix, p.fill,
               p.rows);
}

int chooseCount(double length, double cutoff) {
  int n = static_cast<int>(std::floor(length / cutoff));
  if (n < 1) {
    n = 1;
  }
  return n;
}

int reachOf(double length, int ncell, double cutoff) {
  const double width = length / static_cast<double>(ncell);
  int r = static_cast<int>(std::ceil(cutoff / width - 1e-12));
  if (r < 1) {
    r = 1;
  }
  return r;
}

template <bool Write>
void passFast(const double *xs, const double *ys, const double *zs,
              const int *ids, const int *cellStart, int ncx, int ncy, int ncz,
              double lx, double ly, double lz, double rc2, int *counts,
              int *offsets, int *cols) {
  const int nCells = ncx * ncy * ncz;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nCells >= 64)
#endif
  for (int c = 0; c < nCells; ++c) {
    const int b = cellStart[c];
    const int e = cellStart[c + 1];
    if (b == e) {
      continue;
    }
    const int cz = c % ncz;
    const int cy = (c / ncz) % ncy;
    const int cx = c / (ncy * ncz);
    struct Nbr {
      int begin;
      int end;
      double sx, sy, sz;
    };
    Nbr nbr[27];
    int nn = 0;
    for (int ox = -1; ox <= 1; ++ox) {
      int nx = cx + ox;
      double sx = 0.0;
      if (nx >= ncx) {
        nx -= ncx;
        sx = lx;
      } else if (nx < 0) {
        nx += ncx;
        sx = -lx;
      }
      for (int oy = -1; oy <= 1; ++oy) {
        int ny = cy + oy;
        double sy = 0.0;
        if (ny >= ncy) {
          ny -= ncy;
          sy = ly;
        } else if (ny < 0) {
          ny += ncy;
          sy = -ly;
        }
        const int row = nx * ncy + ny;
        for (int oz = -1; oz <= 1; ++oz) {
          int nz = cz + oz;
          double sz = 0.0;
          if (nz >= ncz) {
            nz -= ncz;
            sz = lz;
          } else if (nz < 0) {
            nz += ncz;
            sz = -lz;
          }
          const int ncid = row * ncz + nz;
          const int jb = cellStart[ncid];
          const int je = cellStart[ncid + 1];
          if (jb == je) {
            continue;
          }
          nbr[nn++] = {jb, je, sx, sy, sz};
        }
      }
    }
    for (int i = b; i < e; ++i) {
      const int ii = ids[i];
      const double xi = xs[i];
      const double yi = ys[i];
      const double zi = zs[i];
      const int base = Write ? offsets[ii] : 0;
      int n = 0;
      for (int t = 0; t < nn; ++t) {
        const Nbr nb = nbr[t];
        for (int j = nb.begin; j < nb.end; ++j) {
          const int jj = ids[j];
          if (jj == ii) {
            continue;
          }
          const double dx = xi - xs[j] - nb.sx;
          const double dy = yi - ys[j] - nb.sy;
          const double dz = zi - zs[j] - nb.sz;
          if (dx * dx + dy * dy + dz * dz < rc2) {
            if (Write) {
              cols[base + n] = jj;
            }
            ++n;
          }
        }
      }
      if (Write) {
        sortRow(cols + base, n);
      } else {
        counts[ii] = n;
      }
    }
  }
}

template <bool Write>
void passSafe(const double *xs, const double *ys, const double *zs,
              const int *ids, const int *cellStart, int ncx, int ncy, int ncz,
              int rx, int ry, int rz, double lx, double ly, double lz,
              double rc2, int *counts, int *offsets, int *cols) {
  const int nCells = ncx * ncy * ncz;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel if (nCells >= 8)
#endif
  {
    std::vector<int> found;
    std::vector<int> stamp(static_cast<std::size_t>(nCells), 0);
    int gen = 0;
    found.reserve(64);
#ifdef SEAMS_HAS_OPENMP
#pragma omp for schedule(static)
#endif
    for (int c = 0; c < nCells; ++c) {
      const int b = cellStart[c];
      const int e = cellStart[c + 1];
      if (b == e) {
        continue;
      }
      const int cz = c % ncz;
      const int cy = (c / ncz) % ncy;
      const int cx = c / (ncy * ncz);
      for (int i = b; i < e; ++i) {
        const int ii = ids[i];
        found.clear();
        ++gen;
        if (gen == 0) {
          std::fill(stamp.begin(), stamp.end(), 0);
          gen = 1;
        }
        for (int ox = -rx; ox <= rx; ++ox) {
          int nx = cx + ox;
          nx %= ncx;
          if (nx < 0) {
            nx += ncx;
          }
          for (int oy = -ry; oy <= ry; ++oy) {
            int ny = cy + oy;
            ny %= ncy;
            if (ny < 0) {
              ny += ncy;
            }
            for (int oz = -rz; oz <= rz; ++oz) {
              int nz = cz + oz;
              nz %= ncz;
              if (nz < 0) {
                nz += ncz;
              }
              const int ncid = (nx * ncy + ny) * ncz + nz;
              if (stamp[static_cast<std::size_t>(ncid)] == gen) {
                continue;
              }
              stamp[static_cast<std::size_t>(ncid)] = gen;
              const int jb = cellStart[ncid];
              const int je = cellStart[ncid + 1];
              for (int j = jb; j < je; ++j) {
                const int jj = ids[j];
                if (jj == ii) {
                  continue;
                }
                const double d2 = minImage2(xs[i] - xs[j], ys[i] - ys[j],
                                           zs[i] - zs[j], lx, ly, lz);
                if (d2 < rc2) {
                  found.push_back(jj);
                }
              }
            }
          }
        }
        std::sort(found.begin(), found.end());
        found.erase(std::unique(found.begin(), found.end()), found.end());
        if (Write) {
          const int base = offsets[ii];
          for (int k = 0; k < static_cast<int>(found.size()); ++k) {
            cols[base + k] = found[static_cast<std::size_t>(k)];
          }
        } else {
          counts[ii] = static_cast<int>(found.size());
        }
      }
    }
  }
}

} // namespace

namespace nneigh {

bool orthoFullCsr(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<int> *subset, double rcutoff, std::vector<int> &offsets,
    std::vector<int> &cols) {
  const int nop = yCloud.nop;
  if (nop < 0 || yCloud.box.size() != 3 || !(rcutoff > 0.0)) {
    return false;
  }
  const double len[3] = {yCloud.box[0], yCloud.box[1], yCloud.box[2]};
  if (!(len[0] > 0.0 && len[1] > 0.0 && len[2] > 0.0)) {
    return false;
  }
  const double origin[3] = {
      yCloud.boxLow.size() > 0 ? yCloud.boxLow[0] : 0.0,
      yCloud.boxLow.size() > 1 ? yCloud.boxLow[1] : 0.0,
      yCloud.boxLow.size() > 2 ? yCloud.boxLow[2] : 0.0};
  const int nSubset = subset == nullptr ? nop : static_cast<int>(subset->size());
  offsets.assign(static_cast<std::size_t>(nop) + 1, 0);
  cols.clear();
  if (nop == 0 || nSubset == 0) {
    return true;
  }
  if (subset != nullptr) {
    for (int k = 0; k < nSubset; ++k) {
      const int id = (*subset)[static_cast<std::size_t>(k)];
      if (id < 0 || id >= nop) {
        return false;
      }
    }
  }

  int ncell[3] = {chooseCount(len[0], rcutoff), chooseCount(len[1], rcutoff),
                  chooseCount(len[2], rcutoff)};
  while (static_cast<long>(ncell[0]) * ncell[1] * ncell[2] > kMaxCells) {
    int *which = &ncell[0];
    if (ncell[1] > *which) {
      which = &ncell[1];
    }
    if (ncell[2] > *which) {
      which = &ncell[2];
    }
    if (*which <= 1) {
      break;
    }
    --(*which);
  }
  const int reach[3] = {reachOf(len[0], ncell[0], rcutoff),
                        reachOf(len[1], ncell[1], rcutoff),
                        reachOf(len[2], ncell[2], rcutoff)};
  const bool fast = reach[0] == 1 && reach[1] == 1 && reach[2] == 1 &&
                    ncell[0] >= 4 && ncell[1] >= 4 && ncell[2] >= 4;
  const int nCells = ncell[0] * ncell[1] * ncell[2];
  const double rc2 = rcutoff * rcutoff;

#ifdef SEAMS_HAS_OPENMP
  omp_set_dynamic(0);
  const int maxThreads = std::max(1, omp_get_max_threads());
#else
  const int maxThreads = 1;
#endif
  Scratch local;
#ifdef SEAMS_HAS_OPENMP
  const bool nested = omp_in_parallel();
#else
  const bool nested = false;
#endif
  Scratch &sc = nested ? local : scratch();
  sc.px.resize(static_cast<std::size_t>(nSubset));
  sc.py.resize(static_cast<std::size_t>(nSubset));
  sc.pz.resize(static_cast<std::size_t>(nSubset));
  sc.xs.resize(static_cast<std::size_t>(nSubset));
  sc.ys.resize(static_cast<std::size_t>(nSubset));
  sc.zs.resize(static_cast<std::size_t>(nSubset));
  sc.cellOf.resize(static_cast<std::size_t>(nSubset));
  sc.ids.resize(static_cast<std::size_t>(nSubset));
  sc.cellStart.assign(static_cast<std::size_t>(nCells) + 1, 0);
  sc.hist.assign(static_cast<std::size_t>(maxThreads) * static_cast<std::size_t>(nCells),
                 0);
  sc.counts.assign(static_cast<std::size_t>(nop), 0);

  const auto tBin = Clock::now();
  int threads = 1;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel if (!nested && nSubset >= 256)
#endif
  {
    const int tid = teamRank();
    const int nt = teamSize();
#ifdef SEAMS_HAS_OPENMP
#pragma omp single
#endif
    {
      threads = nt;
    }
    int *hist = sc.hist.data() + static_cast<std::size_t>(tid) * static_cast<std::size_t>(nCells);
    std::fill(hist, hist + nCells, 0);
#ifdef SEAMS_HAS_OPENMP
#pragma omp for schedule(static)
#endif
    for (int k = 0; k < nSubset; ++k) {
      const int id = subset == nullptr ? k : (*subset)[static_cast<std::size_t>(k)];
      const auto &p = yCloud.pts[static_cast<std::size_t>(id)];
      const double coord[3] = {p.x, p.y, p.z};
      int cxyz[3];
      double wrapped[3];
      for (int d = 0; d < 3; ++d) {
        double f = (coord[d] - origin[d]) / len[d];
        f -= std::floor(f);
        if (!(f < 1.0)) {
          f = 0.0;
        }
        int c = static_cast<int>(f * static_cast<double>(ncell[d]));
        if (c >= ncell[d]) {
          c = ncell[d] - 1;
        }
        if (c < 0) {
          c = 0;
        }
        cxyz[d] = c;
        wrapped[d] = origin[d] + f * len[d];
      }
      const int cid = (cxyz[0] * ncell[1] + cxyz[1]) * ncell[2] + cxyz[2];
      sc.cellOf[static_cast<std::size_t>(k)] = cid;
      sc.px[static_cast<std::size_t>(k)].v = wrapped[0];
      sc.py[static_cast<std::size_t>(k)].v = wrapped[1];
      sc.pz[static_cast<std::size_t>(k)].v = wrapped[2];
      ++hist[cid];
    }
#ifdef SEAMS_HAS_OPENMP
#pragma omp barrier
#pragma omp single
#endif
    {
      sc.cellStart[0] = 0;
      for (int c = 0; c < nCells; ++c) {
        int sum = 0;
        for (int t = 0; t < nt; ++t) {
          sum += sc.hist[static_cast<std::size_t>(t) * static_cast<std::size_t>(nCells) +
                         static_cast<std::size_t>(c)];
        }
        sc.cellStart[static_cast<std::size_t>(c) + 1] =
            sc.cellStart[static_cast<std::size_t>(c)] + sum;
        int acc = sc.cellStart[static_cast<std::size_t>(c)];
        for (int t = 0; t < nt; ++t) {
          int &h = sc.hist[static_cast<std::size_t>(t) * static_cast<std::size_t>(nCells) +
                           static_cast<std::size_t>(c)];
          const int take = h;
          h = acc;
          acc += take;
        }
      }
    }
#ifdef SEAMS_HAS_OPENMP
#pragma omp barrier
#pragma omp for schedule(static)
#endif
    for (int k = 0; k < nSubset; ++k) {
      const int cid = sc.cellOf[static_cast<std::size_t>(k)];
      const int slot = hist[cid]++;
      const int id = subset == nullptr ? k : (*subset)[static_cast<std::size_t>(k)];
      sc.ids[static_cast<std::size_t>(slot)] = id;
      sc.xs[static_cast<std::size_t>(slot)].v = sc.px[static_cast<std::size_t>(k)].v;
      sc.ys[static_cast<std::size_t>(slot)].v = sc.py[static_cast<std::size_t>(k)].v;
      sc.zs[static_cast<std::size_t>(slot)].v = sc.pz[static_cast<std::size_t>(k)].v;
    }
  }
  phase().bin = msSince(tBin);
  phase().threads = threads;
  phase().n = nSubset;

  const double *x = &sc.xs[0].v;
  const double *y = &sc.ys[0].v;
  const double *z = &sc.zs[0].v;

  const auto tCount = Clock::now();
  if (fast) {
    passFast<false>(x, y, z, sc.ids.data(), sc.cellStart.data(), ncell[0], ncell[1],
                    ncell[2], len[0], len[1], len[2], rc2, sc.counts.data(), nullptr,
                    nullptr);
  } else {
    passSafe<false>(x, y, z, sc.ids.data(), sc.cellStart.data(), ncell[0], ncell[1],
                    ncell[2], reach[0], reach[1], reach[2], len[0], len[1], len[2], rc2,
                    sc.counts.data(), nullptr, nullptr);
  }
  phase().count = msSince(tCount);

  const auto tPrefix = Clock::now();
  offsets[0] = 0;
  for (int i = 0; i < nop; ++i) {
    offsets[static_cast<std::size_t>(i) + 1] =
        offsets[static_cast<std::size_t>(i)] + sc.counts[static_cast<std::size_t>(i)];
  }
  const int nnz = offsets[static_cast<std::size_t>(nop)];
  cols.resize(static_cast<std::size_t>(nnz));
  phase().prefix = msSince(tPrefix);
  phase().pairs = nnz;

  const auto tFill = Clock::now();
  if (nnz > 0) {
    if (fast) {
      passFast<true>(x, y, z, sc.ids.data(), sc.cellStart.data(), ncell[0], ncell[1],
                     ncell[2], len[0], len[1], len[2], rc2, nullptr, offsets.data(),
                     cols.data());
    } else {
      passSafe<true>(x, y, z, sc.ids.data(), sc.cellStart.data(), ncell[0], ncell[1],
                     ncell[2], reach[0], reach[1], reach[2], len[0], len[1], len[2], rc2,
                     nullptr, offsets.data(), cols.data());
    }
  }
  phase().fill = msSince(tFill);
  publishProfile();
  return true;
}

std::vector<std::vector<int>> idRowsFromCsr(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<int> &indexToID, const std::vector<int> &offsets,
    const std::vector<int> &cols) {
  const auto t0 = Clock::now();
  const int nop = yCloud.nop;
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(std::max(nop, 0)));
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nop >= 256 && !omp_in_parallel())
#endif
  for (int i = 0; i < nop; ++i) {
    if (static_cast<int>(indexToID.size()) <= i ||
        indexToID[static_cast<std::size_t>(i)] == -1) {
      continue;
    }
    const int a = offsets[static_cast<std::size_t>(i)];
    const int b = offsets[static_cast<std::size_t>(i) + 1];
    auto &row = rows[static_cast<std::size_t>(i)];
    row.reserve(static_cast<std::size_t>(1 + b - a));
    row.push_back(indexToID[static_cast<std::size_t>(i)]);
    for (int k = a; k < b; ++k) {
      const int j = cols[static_cast<std::size_t>(k)];
      if (j >= 0 && j < static_cast<int>(indexToID.size()) &&
          indexToID[static_cast<std::size_t>(j)] != -1) {
        row.push_back(indexToID[static_cast<std::size_t>(j)]);
      }
    }
  }
  phase().rows = msSince(t0);
  if (profileOn()) {
    publishProfile();
  }
  return rows;
}

std::vector<std::vector<int>> indexRowsFromCsr(int nop, const std::vector<int> &offsets,
                                              const std::vector<int> &cols) {
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(std::max(nop, 0)));
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nop >= 256 && !omp_in_parallel())
#endif
  for (int i = 0; i < nop; ++i) {
    const int a = offsets[static_cast<std::size_t>(i)];
    const int b = offsets[static_cast<std::size_t>(i) + 1];
    auto &row = rows[static_cast<std::size_t>(i)];
    row.resize(static_cast<std::size_t>(1 + b - a));
    row[0] = i;
    for (int k = a; k < b; ++k) {
      row[static_cast<std::size_t>(1 + k - a)] = cols[static_cast<std::size_t>(k)];
    }
  }
  return rows;
}

std::vector<std::vector<int>> halfIdRowsFromCsr(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<int> &indexToID, const std::vector<int> &offsets,
    const std::vector<int> &cols) {
  const int nop = yCloud.nop;
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(std::max(nop, 0)));
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nop >= 256 && !omp_in_parallel())
#endif
  for (int i = 0; i < nop; ++i) {
    if (static_cast<int>(indexToID.size()) <= i ||
        indexToID[static_cast<std::size_t>(i)] == -1) {
      continue;
    }
    const int a = offsets[static_cast<std::size_t>(i)];
    const int b = offsets[static_cast<std::size_t>(i) + 1];
    auto &row = rows[static_cast<std::size_t>(i)];
    row.push_back(indexToID[static_cast<std::size_t>(i)]);
    for (int k = a; k < b; ++k) {
      const int j = cols[static_cast<std::size_t>(k)];
      if (j <= i) {
        continue;
      }
      if (j >= 0 && j < static_cast<int>(indexToID.size()) &&
          indexToID[static_cast<std::size_t>(j)] != -1) {
        row.push_back(indexToID[static_cast<std::size_t>(j)]);
      }
    }
  }
  return rows;
}

} // namespace nneigh
