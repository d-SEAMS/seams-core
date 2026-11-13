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
#include <mutex>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

#ifdef SEAMS_HAS_LINKCELL
#include <linkcell.h>
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
  // Neighbour ids captured during the count walk, stride kStashCap.
  // The fill then copies this instead of walking the cells again.
  std::vector<int> stash;
  int overflowed = 0;
};

constexpr int kStashCap = 48;

Scratch &scratch() {
  // One buffer per calling thread. A trial pool outside OpenMP used to
  // share this scratch and tear the neighbour list.
  thread_local Scratch s;
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
  thread_local Phase p;
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

// shx/shy/shz hold the cartesian image of fractional shift (ix,iy,iz),
// indexed by (ix+1)*9 + (iy+1)*3 + (iz+1), ix,iy,iz in {-1,0,1}.
void fillOrthoShift(double lx, double ly, double lz, double *shx, double *shy,
                   double *shz) {
  for (int ix = -1; ix <= 1; ++ix) {
    for (int iy = -1; iy <= 1; ++iy) {
      for (int iz = -1; iz <= 1; ++iz) {
        const int k = (ix + 1) * 9 + (iy + 1) * 3 + (iz + 1);
        shx[k] = lx * static_cast<double>(ix);
        shy[k] = ly * static_cast<double>(iy);
        shz[k] = lz * static_cast<double>(iz);
      }
    }
  }
}

template <bool Write, bool Half>
void passFast(const double *xs, const double *ys, const double *zs,
              const int *ids, const int *cellStart, int ncx, int ncy, int ncz,
              const double *shx, const double *shy, const double *shz, double rc2,
              int *counts, int *offsets, int *cols) {
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
      int ix = 0;
      if (nx >= ncx) {
        nx -= ncx;
        ix = 1;
      } else if (nx < 0) {
        nx += ncx;
        ix = -1;
      }
      for (int oy = -1; oy <= 1; ++oy) {
        int ny = cy + oy;
        int iy = 0;
        if (ny >= ncy) {
          ny -= ncy;
          iy = 1;
        } else if (ny < 0) {
          ny += ncy;
          iy = -1;
        }
        const int row = nx * ncy + ny;
        for (int oz = -1; oz <= 1; ++oz) {
          int nz = cz + oz;
          int iz = 0;
          if (nz >= ncz) {
            nz -= ncz;
            iz = 1;
          } else if (nz < 0) {
            nz += ncz;
            iz = -1;
          }
          const int ncid = row * ncz + nz;
          const int jb = cellStart[ncid];
          const int je = cellStart[ncid + 1];
          if (jb == je) {
            continue;
          }
          const int k = (ix + 1) * 9 + (iy + 1) * 3 + (iz + 1);
          nbr[nn++] = {jb, je, shx[k], shy[k], shz[k]};
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
#pragma GCC unroll 1
        for (int j = nb.begin; j < nb.end; ++j) {
          // Slot identity: each atom is stored once, so j == i is the self pair.
          // The full-list count then does not touch the id stream.
          if (j == i || (Half && ids[j] < ii)) {
            continue;
          }
          const double dx = xi - xs[j] - nb.sx;
          const double dy = yi - ys[j] - nb.sy;
          const double dz = zi - zs[j] - nb.sz;
          if (dx * dx + dy * dy + dz * dz < rc2) {
            if (Write) {
              cols[base + n] = ids[j];
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

// One z layer of the half stencil. Layers two apart do not write the
// same atom: a layer records its own atoms and the atoms in the next
// layer, and the next layer is the other parity. Even layers run
// together, then odd layers.
void halfZLayer(int cz, const double *xs, const double *ys, const double *zs,
                const int *ids, const int *cellStart, int ncx, int ncy, int ncz,
                const double *shx, const double *shy, const double *shz, double rc2,
                int *counts, int *stash, int *overflowed) {
  for (int cy = 0; cy < ncy; ++cy) {
    for (int cx = 0; cx < ncx; ++cx) {
      const int c = (cx * ncy + cy) * ncz + cz;
      const int b = cellStart[c];
      const int e = cellStart[c + 1];
      if (b == e) {
        continue;
      }
      // Interior cells have no periodic image. The thirteen forward
      // neighbours are plain index steps, and the shift is zero.
      if (cx > 0 && cx + 1 < ncx && cy > 0 && cy + 1 < ncy && cz > 0 && cz + 1 < ncz) {
        struct Step {
          int begin;
          int end;
        };
        Step nbr[14];
        int nn = 0;
        nbr[nn++] = {b, e};
        const int delta[13][3] = {
            {1, -1, -1}, {1, -1, 0}, {1, -1, 1}, {1, 0, -1}, {1, 0, 0},
            {1, 0, 1},   {1, 1, -1}, {1, 1, 0},  {1, 1, 1},  {0, 1, -1},
            {0, 1, 0},   {0, 1, 1},   {0, 0, 1},
        };
        for (const auto &d : delta) {
          const int ncid = ((cx + d[0]) * ncy + (cy + d[1])) * ncz + (cz + d[2]);
          const int jb = cellStart[ncid];
          const int je = cellStart[ncid + 1];
          if (jb != je) {
            nbr[nn++] = {jb, je};
          }
        }
        for (int i = b; i < e; ++i) {
          const int ii = ids[i];
          const double xi = xs[i];
          const double yi = ys[i];
          const double zi = zs[i];
          int *rowI = stash + static_cast<std::size_t>(ii) * static_cast<std::size_t>(kStashCap);
          for (int t = 0; t < nn; ++t) {
            const int j0 = (t == 0) ? i + 1 : nbr[t].begin;
            for (int j = j0; j < nbr[t].end; ++j) {
              const double dx = xi - xs[j];
              const double dy = yi - ys[j];
              const double dz = zi - zs[j];
              if (!(dx * dx + dy * dy + dz * dz < rc2)) {
                continue;
              }
              const int jj = ids[j];
              int &ni = counts[ii];
              int &nj = counts[jj];
              if (ni < kStashCap) {
                rowI[ni] = jj;
              } else {
#ifdef SEAMS_HAS_OPENMP
#pragma omp atomic write
#endif
                *overflowed = 1;
              }
              ++ni;
              int *rowJ = stash + static_cast<std::size_t>(jj) * static_cast<std::size_t>(kStashCap);
              if (nj < kStashCap) {
                rowJ[nj] = ii;
              } else {
#ifdef SEAMS_HAS_OPENMP
#pragma omp atomic write
#endif
                *overflowed = 1;
              }
              ++nj;
            }
          }
        }
        continue;
      }
      struct Nbr {
        int begin;
        int end;
        double sx, sy, sz;
      };
      Nbr nbr[14];
      int nn = 0;
      nbr[nn++] = {b, e, 0.0, 0.0, 0.0};
      for (int ox = -1; ox <= 1; ++ox) {
        for (int oy = -1; oy <= 1; ++oy) {
          for (int oz = -1; oz <= 1; ++oz) {
            if (!(ox > 0 || (ox == 0 && oy > 0) || (ox == 0 && oy == 0 && oz > 0))) {
              continue;
            }
            int nx = cx + ox;
            int ix = 0;
            if (nx >= ncx) {
              nx -= ncx;
              ix = 1;
            } else if (nx < 0) {
              nx += ncx;
              ix = -1;
            }
            int ny = cy + oy;
            int iy = 0;
            if (ny >= ncy) {
              ny -= ncy;
              iy = 1;
            } else if (ny < 0) {
              ny += ncy;
              iy = -1;
            }
            int nz = cz + oz;
            int iz = 0;
            if (nz >= ncz) {
              nz -= ncz;
              iz = 1;
            } else if (nz < 0) {
              nz += ncz;
              iz = -1;
            }
            const int ncid = (nx * ncy + ny) * ncz + nz;
            const int jb = cellStart[ncid];
            const int je = cellStart[ncid + 1];
            if (jb == je) {
              continue;
            }
            const int k = (ix + 1) * 9 + (iy + 1) * 3 + (iz + 1);
            nbr[nn++] = {jb, je, shx[k], shy[k], shz[k]};
          }
        }
      }
      for (int i = b; i < e; ++i) {
        const int ii = ids[i];
        const double xi = xs[i];
        const double yi = ys[i];
        const double zi = zs[i];
        int *rowI = stash + static_cast<std::size_t>(ii) * static_cast<std::size_t>(kStashCap);
        for (int t = 0; t < nn; ++t) {
          const Nbr nb = nbr[t];
          const int j0 = (t == 0) ? i + 1 : nb.begin;
          for (int j = j0; j < nb.end; ++j) {
            const double dx = xi - xs[j] - nb.sx;
            const double dy = yi - ys[j] - nb.sy;
            const double dz = zi - zs[j] - nb.sz;
            if (!(dx * dx + dy * dy + dz * dz < rc2)) {
              continue;
            }
            const int jj = ids[j];
            int &ni = counts[ii];
            int &nj = counts[jj];
            if (ni < kStashCap) {
              rowI[ni] = jj;
            } else {
#ifdef SEAMS_HAS_OPENMP
#pragma omp atomic write
#endif
              *overflowed = 1;
            }
            ++ni;
            int *rowJ = stash + static_cast<std::size_t>(jj) * static_cast<std::size_t>(kStashCap);
            if (nj < kStashCap) {
              rowJ[nj] = ii;
            } else {
#ifdef SEAMS_HAS_OPENMP
#pragma omp atomic write
#endif
              *overflowed = 1;
            }
            ++nj;
          }
        }
      }
    }
  }
}

void passFastStashHalf(const double *xs, const double *ys, const double *zs,
                       const int *ids, const int *cellStart, int ncx, int ncy, int ncz,
                       const double *shx, const double *shy, const double *shz, double rc2,
                       int *counts, int *stash, int *overflowed);

// Count pass that also keeps the neighbour ids. A shell longer than
// kStashCap sets *overflowed and the caller repeats the geometry fill.
// Several threads take the half stencil a z layer at a time. The full
// stencil visited every pair twice and did not beat one thread.
void passFastStash(const double *xs, const double *ys, const double *zs,
                   const int *ids, const int *cellStart, int ncx, int ncy, int ncz,
                   const double *shx, const double *shy, const double *shz, double rc2,
                   int *counts, int *stash, int *overflowed) {
  *overflowed = 0;
#ifdef SEAMS_HAS_OPENMP
  if (ncz >= 2 && omp_get_max_threads() > 1) {
    // A layer writes its own atoms and the next layer's. Layers of one
    // parity therefore do not share an atom. When the length is odd the
    // last layer's +z wraps onto layer 0, which is the same parity, so
    // that layer runs alone after the others have stopped.
    const int lastWraps = (ncz % 2) == 1;
    const int nEven = ncz / 2;
#pragma omp parallel for schedule(static)
    for (int k = 0; k < nEven; ++k) {
      halfZLayer(k * 2, xs, ys, zs, ids, cellStart, ncx, ncy, ncz, shx, shy, shz, rc2,
                 counts, stash, overflowed);
    }
    if (lastWraps) {
      halfZLayer(ncz - 1, xs, ys, zs, ids, cellStart, ncx, ncy, ncz, shx, shy, shz, rc2,
                 counts, stash, overflowed);
    }
    const int nOdd = ncz / 2;
#pragma omp parallel for schedule(static)
    for (int k = 0; k < nOdd; ++k) {
      halfZLayer(k * 2 + 1, xs, ys, zs, ids, cellStart, ncx, ncy, ncz, shx, shy, shz, rc2,
                 counts, stash, overflowed);
    }
    return;
  }
#endif
  passFastStashHalf(xs, ys, zs, ids, cellStart, ncx, ncy, ncz, shx, shy, shz, rc2, counts,
                    stash, overflowed);
}

// One thread records each pair from the lower packed slot and writes both
// ends. The full stencil visits every pair twice; this visit is once.
void passFastStashHalf(const double *xs, const double *ys, const double *zs,
                       const int *ids, const int *cellStart, int ncx, int ncy, int ncz,
                       const double *shx, const double *shy, const double *shz, double rc2,
                       int *counts, int *stash, int *overflowed) {
  const int nCells = ncx * ncy * ncz;
  *overflowed = 0;
  for (int c = 0; c < nCells; ++c) {
    const int b = cellStart[c];
    const int e = cellStart[c + 1];
    if (b == e) {
      continue;
    }
    const int cz = c % ncz;
    const int cy = (c / ncz) % ncy;
    const int cx = c / (ncy * ncz);
    // An interior cell has no periodic image. The forward neighbours are
    // index steps and the shift is zero, so the wrap tests stay on the
    // boundary layer.
    if (cx > 0 && cx + 1 < ncx && cy > 0 && cy + 1 < ncy && cz > 0 && cz + 1 < ncz) {
      int begin[14];
      int end[14];
      int nn = 0;
      begin[nn] = b;
      end[nn] = e;
      ++nn;
      const int delta[13][3] = {
          {1, -1, -1}, {1, -1, 0}, {1, -1, 1}, {1, 0, -1}, {1, 0, 0},
          {1, 0, 1},   {1, 1, -1}, {1, 1, 0},  {1, 1, 1},  {0, 1, -1},
          {0, 1, 0},   {0, 1, 1},   {0, 0, 1},
      };
      for (const auto &d : delta) {
        const int ncid = ((cx + d[0]) * ncy + (cy + d[1])) * ncz + (cz + d[2]);
        const int jb = cellStart[ncid];
        const int je = cellStart[ncid + 1];
        if (jb != je) {
          begin[nn] = jb;
          end[nn] = je;
          ++nn;
        }
      }
      for (int i = b; i < e; ++i) {
        const int ii = ids[i];
        const double xi = xs[i];
        const double yi = ys[i];
        const double zi = zs[i];
        int *rowI = stash + static_cast<std::size_t>(ii) * static_cast<std::size_t>(kStashCap);
        for (int t = 0; t < nn; ++t) {
          const int j0 = (t == 0) ? i + 1 : begin[t];
          for (int j = j0; j < end[t]; ++j) {
            const double dx = xi - xs[j];
            const double dy = yi - ys[j];
            const double dz = zi - zs[j];
            if (!(dx * dx + dy * dy + dz * dz < rc2)) {
              continue;
            }
            const int jj = ids[j];
            int &ni = counts[ii];
            int &nj = counts[jj];
            if (ni < kStashCap) {
              rowI[ni] = jj;
            } else {
              *overflowed = 1;
            }
            ++ni;
            int *rowJ = stash + static_cast<std::size_t>(jj) * static_cast<std::size_t>(kStashCap);
            if (nj < kStashCap) {
              rowJ[nj] = ii;
            } else {
              *overflowed = 1;
            }
            ++nj;
          }
        }
      }
      continue;
    }
    struct Nbr {
      int begin;
      int end;
      double sx, sy, sz;
    };
    Nbr nbr[14];
    int nn = 0;
    // Same cell, then the 13 cells on the positive side of the stencil.
    // The opposite cell records the pair, so each pair is seen once.
    nbr[nn++] = {b, e, 0.0, 0.0, 0.0};
    for (int ox = -1; ox <= 1; ++ox) {
      for (int oy = -1; oy <= 1; ++oy) {
        for (int oz = -1; oz <= 1; ++oz) {
          if (!(ox > 0 || (ox == 0 && oy > 0) || (ox == 0 && oy == 0 && oz > 0))) {
            continue;
          }
          int nx = cx + ox;
          int ix = 0;
          if (nx >= ncx) {
            nx -= ncx;
            ix = 1;
          } else if (nx < 0) {
            nx += ncx;
            ix = -1;
          }
          int ny = cy + oy;
          int iy = 0;
          if (ny >= ncy) {
            ny -= ncy;
            iy = 1;
          } else if (ny < 0) {
            ny += ncy;
            iy = -1;
          }
          int nz = cz + oz;
          int iz = 0;
          if (nz >= ncz) {
            nz -= ncz;
            iz = 1;
          } else if (nz < 0) {
            nz += ncz;
            iz = -1;
          }
          const int ncid = (nx * ncy + ny) * ncz + nz;
          const int jb = cellStart[ncid];
          const int je = cellStart[ncid + 1];
          if (jb == je) {
            continue;
          }
          const int k = (ix + 1) * 9 + (iy + 1) * 3 + (iz + 1);
          nbr[nn++] = {jb, je, shx[k], shy[k], shz[k]};
        }
      }
    }
    for (int i = b; i < e; ++i) {
      const int ii = ids[i];
      const double xi = xs[i];
      const double yi = ys[i];
      const double zi = zs[i];
      int *rowI = stash + static_cast<std::size_t>(ii) * static_cast<std::size_t>(kStashCap);
      for (int t = 0; t < nn; ++t) {
        const Nbr nb = nbr[t];
        const int j0 = (t == 0) ? i + 1 : nb.begin;
        for (int j = j0; j < nb.end; ++j) {
          const double dx = xi - xs[j] - nb.sx;
          const double dy = yi - ys[j] - nb.sy;
          const double dz = zi - zs[j] - nb.sz;
          if (!(dx * dx + dy * dy + dz * dz < rc2)) {
            continue;
          }
          const int jj = ids[j];
          int &ni = counts[ii];
          int &nj = counts[jj];
          if (ni < kStashCap) {
            rowI[ni] = jj;
          } else {
            *overflowed = 1;
          }
          ++ni;
          int *rowJ = stash + static_cast<std::size_t>(jj) * static_cast<std::size_t>(kStashCap);
          if (nj < kStashCap) {
            rowJ[nj] = ii;
          } else {
            *overflowed = 1;
          }
          ++nj;
        }
      }
    }
  }
}

void scatterStash(const int *stash, const int *counts, const int *offsets, int *cols,
                  int nop) {
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nop >= 256)
#endif
  for (int i = 0; i < nop; ++i) {
    const int n = counts[i];
    const int *src = stash + static_cast<std::size_t>(i) * static_cast<std::size_t>(kStashCap);
    int *dst = cols + offsets[i];
    for (int k = 0; k < n; ++k) {
      dst[k] = src[k];
    }
    sortRow(dst, n);
  }
}

template <bool Write>
void passSafe(const double *xs, const double *ys, const double *zs,
              const int *ids, const int *cellStart, int ncx, int ncy, int ncz,
              int rx, int ry, int rz, double lx, double ly, double lz,
              double rc2, int *counts, int *offsets, int *cols, bool half) {
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
                if (jj == ii || (half && jj < ii)) {
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

struct Tilt {
  double lx, ly, lz, xy, xz, yz, ox, oy, oz;
};

double tilt2(double dsx, double dsy, double dsz, const Tilt &h) {
  const double dx = h.lx * dsx + h.xy * dsy + h.xz * dsz;
  const double dy = h.ly * dsy + h.yz * dsz;
  const double dz = h.lz * dsz;
  return dx * dx + dy * dy + dz * dz;
}

void fracOf(const Tilt &h, double x, double y, double z, double s[3]) {
  s[2] = (z - h.oz) / h.lz;
  s[1] = (y - h.oy - h.yz * s[2]) / h.ly;
  s[0] = (x - h.ox - h.xy * s[1] - h.xz * s[2]) / h.lx;
  for (int d = 0; d < 3; ++d) {
    s[d] -= std::floor(s[d]);
    if (!(s[d] < 1.0)) {
      s[d] = 0.0;
    }
  }
}

template <bool Write>
void passTilt(const double *xs, const double *ys, const double *zs, const int *ids,
              const int *cellStart, int ncx, int ncy, int ncz, int rx, int ry, int rz,
              const Tilt &h, double rc2, bool half, bool shifted, int *counts,
              int *offsets, int *cols) {
  const int nCells = ncx * ncy * ncz;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel if (nCells >= 8)
#endif
  {
    std::vector<int> found;
    std::vector<int> stamp(static_cast<std::size_t>(nCells), 0);
    int gen = 0;
    found.reserve(32);
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
          double ssx = 0.0;
          while (nx >= ncx) {
            nx -= ncx;
            ssx += 1.0;
          }
          while (nx < 0) {
            nx += ncx;
            ssx -= 1.0;
          }
          for (int oy = -ry; oy <= ry; ++oy) {
            int ny = cy + oy;
            double ssy = 0.0;
            while (ny >= ncy) {
              ny -= ncy;
              ssy += 1.0;
            }
            while (ny < 0) {
              ny += ncy;
              ssy -= 1.0;
            }
            for (int oz = -rz; oz <= rz; ++oz) {
              int nz = cz + oz;
              double ssz = 0.0;
              while (nz >= ncz) {
                nz -= ncz;
                ssz += 1.0;
              }
              while (nz < 0) {
                nz += ncz;
                ssz -= 1.0;
              }
              const int ncid = (nx * ncy + ny) * ncz + nz;
              if (stamp[static_cast<std::size_t>(ncid)] == gen) {
                continue;
              }
              stamp[static_cast<std::size_t>(ncid)] = gen;
              for (int j = cellStart[ncid]; j < cellStart[ncid + 1]; ++j) {
                const int jj = ids[j];
                if (jj == ii || (half && jj < ii)) {
                  continue;
                }
                double dsx, dsy, dsz;
                if (shifted) {
                  dsx = xs[i] - xs[j] - ssx;
                  dsy = ys[i] - ys[j] - ssy;
                  dsz = zs[i] - zs[j] - ssz;
                } else {
                  dsx = xs[i] - xs[j];
                  dsy = ys[i] - ys[j];
                  dsz = zs[i] - zs[j];
                  dsx -= std::floor(dsx + 0.5);
                  dsy -= std::floor(dsy + 0.5);
                  dsz -= std::floor(dsz + 0.5);
                }
                if (tilt2(dsx, dsy, dsz, h) < rc2) {
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

bool tiltFullCsr(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                 const std::vector<int> *subset, double rcutoff,
                 std::vector<int> &offsets, std::vector<int> &cols, bool half) {
  const gen::FracBox box = gen::makeFracBox(yCloud);
  if (!box.ok || !box.triclinic) {
    return false;
  }
  const int nop = yCloud.nop;
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
  const Tilt h{box.lx, box.ly, box.lz, box.xy, box.xz, box.yz, box.ox, box.oy, box.oz};
  const double width[3] = {box.wx, box.wy, box.wz};
  int ncell[3];
  for (int d = 0; d < 3; ++d) {
    ncell[d] = width[d] > 0.0 ? static_cast<int>(std::floor(width[d] / rcutoff)) : 1;
    if (ncell[d] < 1) {
      ncell[d] = 1;
    }
  }
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
  int reach[3];
  for (int d = 0; d < 3; ++d) {
    const double cellW = width[d] / static_cast<double>(ncell[d]);
    reach[d] = cellW > 0.0 ? static_cast<int>(std::ceil(rcutoff / cellW - 1e-12)) : 1;
    if (reach[d] < 1) {
      reach[d] = 1;
    }
  }
  // Two cells along an edge are enough for the one-cell stencil. The
  // other image of a neighbour then sits at least half an edge away, and
  // half an edge is at least the cutoff, so the strict test keeps one
  // image. One cell would list several images of the same atom.
  const bool shifted = reach[0] == 1 && reach[1] == 1 && reach[2] == 1 && ncell[0] >= 2 &&
                       ncell[1] >= 2 && ncell[2] >= 2;
  const int nCells = ncell[0] * ncell[1] * ncell[2];
  const double rc2 = rcutoff * rcutoff;
#ifdef SEAMS_HAS_OPENMP
  omp_set_dynamic(0);
  const int maxThreads = std::max(1, omp_get_max_threads());
  const bool nested = omp_in_parallel();
#else
  const int maxThreads = 1;
  const bool nested = false;
#endif
  Scratch local;
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
  sc.hist.assign(static_cast<std::size_t>(maxThreads) * static_cast<std::size_t>(nCells), 0);
  sc.counts.assign(static_cast<std::size_t>(nop), 0);
  const auto tBin = Clock::now();
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel if (!nested && nSubset >= 256)
#endif
  {
    const int tid = teamRank();
    const int nt = teamSize();
    int *hist = sc.hist.data() + static_cast<std::size_t>(tid) * static_cast<std::size_t>(nCells);
    std::fill(hist, hist + nCells, 0);
#ifdef SEAMS_HAS_OPENMP
#pragma omp for schedule(static)
#endif
    for (int k = 0; k < nSubset; ++k) {
      const int id = subset == nullptr ? k : (*subset)[static_cast<std::size_t>(k)];
      const auto &p = yCloud.pts[static_cast<std::size_t>(id)];
      double s[3];
      fracOf(h, p.x, p.y, p.z, s);
      int cxyz[3];
      for (int d = 0; d < 3; ++d) {
        int c = static_cast<int>(s[d] * static_cast<double>(ncell[d]));
        if (c >= ncell[d]) {
          c = ncell[d] - 1;
        }
        if (c < 0) {
          c = 0;
        }
        cxyz[d] = c;
      }
      const int cid = (cxyz[0] * ncell[1] + cxyz[1]) * ncell[2] + cxyz[2];
      sc.cellOf[static_cast<std::size_t>(k)] = cid;
      sc.px[static_cast<std::size_t>(k)].v = s[0];
      sc.py[static_cast<std::size_t>(k)].v = s[1];
      sc.pz[static_cast<std::size_t>(k)].v = s[2];
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
          int &slot = sc.hist[static_cast<std::size_t>(t) * static_cast<std::size_t>(nCells) +
                              static_cast<std::size_t>(c)];
          const int take = slot;
          slot = acc;
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
      const double s0 = sc.px[static_cast<std::size_t>(k)].v;
      const double s1 = sc.py[static_cast<std::size_t>(k)].v;
      const double s2 = sc.pz[static_cast<std::size_t>(k)].v;
      if (shifted) {
        sc.xs[static_cast<std::size_t>(slot)].v = h.ox + h.lx * s0 + h.xy * s1 + h.xz * s2;
        sc.ys[static_cast<std::size_t>(slot)].v = h.oy + h.ly * s1 + h.yz * s2;
        sc.zs[static_cast<std::size_t>(slot)].v = h.oz + h.lz * s2;
      } else {
        // The wide-reach walk still takes the fractional minimum image.
        sc.xs[static_cast<std::size_t>(slot)].v = s0;
        sc.ys[static_cast<std::size_t>(slot)].v = s1;
        sc.zs[static_cast<std::size_t>(slot)].v = s2;
      }
    }
  }
  phase().bin = msSince(tBin);
  phase().n = nSubset;
  const double *x = &sc.xs[0].v;
  const double *y = &sc.ys[0].v;
  const double *z = &sc.zs[0].v;
  double shx[27], shy[27], shz[27];
  for (int ix = -1; ix <= 1; ++ix) {
    for (int iy = -1; iy <= 1; ++iy) {
      for (int iz = -1; iz <= 1; ++iz) {
        const int k = (ix + 1) * 9 + (iy + 1) * 3 + (iz + 1);
        shx[k] = h.lx * static_cast<double>(ix) + h.xy * static_cast<double>(iy) +
                 h.xz * static_cast<double>(iz);
        shy[k] = h.ly * static_cast<double>(iy) + h.yz * static_cast<double>(iz);
        shz[k] = h.lz * static_cast<double>(iz);
      }
    }
  }
  const auto tCount = Clock::now();
  if (shifted) {
    if (half) {
      passFast<false, true>(x, y, z, sc.ids.data(), sc.cellStart.data(), ncell[0], ncell[1],
                            ncell[2], shx, shy, shz, rc2, sc.counts.data(), nullptr, nullptr);
    } else {
      passFast<false, false>(x, y, z, sc.ids.data(), sc.cellStart.data(), ncell[0], ncell[1],
                             ncell[2], shx, shy, shz, rc2, sc.counts.data(), nullptr, nullptr);
    }
  } else {
    passTilt<false>(x, y, z, sc.ids.data(), sc.cellStart.data(), ncell[0], ncell[1], ncell[2],
                    reach[0], reach[1], reach[2], h, rc2, half, false, sc.counts.data(),
                    nullptr, nullptr);
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
  const auto tFill = Clock::now();
  if (nnz > 0) {
    if (shifted) {
      if (half) {
        passFast<true, true>(x, y, z, sc.ids.data(), sc.cellStart.data(), ncell[0], ncell[1],
                             ncell[2], shx, shy, shz, rc2, nullptr, offsets.data(), cols.data());
      } else {
        passFast<true, false>(x, y, z, sc.ids.data(), sc.cellStart.data(), ncell[0], ncell[1],
                              ncell[2], shx, shy, shz, rc2, nullptr, offsets.data(), cols.data());
      }
    } else {
      passTilt<true>(x, y, z, sc.ids.data(), sc.cellStart.data(), ncell[0], ncell[1], ncell[2],
                     reach[0], reach[1], reach[2], h, rc2, half, false, nullptr, offsets.data(),
                     cols.data());
    }
  }
  phase().fill = msSince(tFill);
  phase().pairs = nnz;
  if (profileOn()) {
    std::fprintf(stderr, "tilt shifted=%d ncell=%d %d %d reach=%d %d %d\n", int(shifted),
                 ncell[0], ncell[1], ncell[2], reach[0], reach[1], reach[2]);
    publishProfile();
  }
  return true;
}

} // namespace

namespace nneigh {

struct PairDir {
  int a;
  int b;
};

#ifdef SEAMS_HAS_LINKCELL
int ourToLc(int cid, int ncx, int ncy, int ncz) {
  const int nyz = ncy * ncz;
  const int cx = cid / nyz;
  const int cy = (cid / ncz) % ncy;
  const int cz = cid % ncz;
  return (cz * ncy + cy) * ncx + cx;
}

// Each thread owns one part of the linkcell bisection and keeps a private
// copy of the atoms in the cells it owns or paints. The search never walks
// another thread's interior.
bool ownedBlockPairs(const double *x, const double *y, const double *z, const int *ids,
                     const int *cellStart, const int *cellOfSubset, const std::vector<int> *subset,
                     int ncx, int ncy, int ncz, int nop, int nSubset, double cutoff, double rc2,
                     const double len[3], int team, std::vector<std::vector<PairDir>> &bags,
                     std::vector<int> &deg, std::vector<std::vector<int>> *directRows,
                     const int *indexToID) {
  const int nCells = ncx * ncy * ncz;
  if (nCells <= 0 || team < 2 || nSubset <= 0) {
    return false;
  }
  const int nxyz[3] = {ncx, ncy, ncz};
  const double bw[3] = {len[0] / static_cast<double>(ncx), len[1] / static_cast<double>(ncy),
                        len[2] / static_cast<double>(ncz)};
  std::uint32_t hw[3] = {1, 1, 1};
  if (lc_halo_width(cutoff, 1, bw, hw) != 0) {
    return false;
  }
  struct CachedCut {
    int ncx = 0;
    int ncy = 0;
    int ncz = 0;
    int team = 0;
    double lx = 0;
    double ly = 0;
    double lz = 0;
    double cutoff = 0;
    std::vector<int> assign;
  };
  static CachedCut cached;
  static std::mutex cutMu;
  std::vector<int> assign(static_cast<std::size_t>(nCells), 0);
  {
    std::lock_guard<std::mutex> lock(cutMu);
    const bool same = cached.ncx == ncx && cached.ncy == ncy && cached.ncz == ncz &&
                      cached.team == team && cached.lx == len[0] && cached.ly == len[1] &&
                      cached.lz == len[2] && cached.cutoff == cutoff &&
                      static_cast<int>(cached.assign.size()) == nCells;
    if (same) {
      assign = cached.assign;
    } else {
      std::vector<std::uint64_t> weight(static_cast<std::size_t>(nCells), 0);
      for (int c = 0; c < nCells; ++c) {
        const int n = cellStart[c + 1] - cellStart[c];
        if (n > 0) {
          weight[static_cast<std::size_t>(ourToLc(c, ncx, ncy, ncz))] =
              static_cast<std::uint64_t>(n);
        }
      }
      const auto cut0 = Clock::now();
      if (lc_partition_weights(nxyz, weight.data(), static_cast<std::size_t>(nCells),
                               static_cast<std::size_t>(team), 2, 0.0, assign.data()) != 0) {
        return false;
      }
      if (profileOn()) {
        std::fprintf(stderr, "rcb_ms %.3f cells %d\n", msSince(cut0), nCells);
      }
      cached.ncx = ncx;
      cached.ncy = ncy;
      cached.ncz = ncz;
      cached.team = team;
      cached.lx = len[0];
      cached.ly = len[1];
      cached.lz = len[2];
      cached.cutoff = cutoff;
      cached.assign = assign;
    }
  }
  (void)cellOfSubset;
  (void)subset;
  bags.assign(static_cast<std::size_t>(team), {});
  if (directRows == nullptr) {
    const std::size_t guess =
        static_cast<std::size_t>(nSubset) * 8 / static_cast<std::size_t>(team) + 64;
    for (auto &bag : bags) {
      bag.reserve(guess);
    }
  }
  if (directRows == nullptr) {
    deg.assign(static_cast<std::size_t>(team) * static_cast<std::size_t>(nop), 0);
  } else {
    deg.clear();
    if (static_cast<int>(directRows->size()) != nop) {
      directRows->assign(static_cast<std::size_t>(nop), {});
      for (auto &row : *directRows) {
        row.reserve(8);
      }
    }
  }
  const int rx = static_cast<int>(hw[0]);
  const int ry = static_cast<int>(hw[1]);
  const int rz = static_cast<int>(hw[2]);
  const bool near = rx == 1 && ry == 1 && rz == 1;
  double shx[27], shy[27], shz[27];
  fillOrthoShift(len[0], len[1], len[2], shx, shy, shz);
  int badTeam = 0;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel num_threads(team)
#endif
  {
    int tid = 0;
    int nt = 1;
#ifdef SEAMS_HAS_OPENMP
    tid = omp_get_thread_num();
    nt = omp_get_num_threads();
#endif
    if (nt != team) {
#ifdef SEAMS_HAS_OPENMP
#pragma omp atomic write
#endif
      badTeam = 1;
    } else {
      const int nyz = ncy * ncz;
      const auto block0 = Clock::now();
      std::vector<int> owned;
      owned.reserve(static_cast<std::size_t>(nCells) / static_cast<std::size_t>(team) + 16);
      for (int c = 0; c < nCells; ++c) {
        if (cellStart[c] == cellStart[c + 1]) {
          continue;
        }
        if (assign[static_cast<std::size_t>(ourToLc(c, ncx, ncy, ncz))] == tid) {
          owned.push_back(c);
        }
      }
      // 2 = owned, 1 = halo painted while walking the owned stencil, 0 = outside.
      std::vector<char> mark(static_cast<std::size_t>(nCells), 0);
      for (int c : owned) {
        mark[static_cast<std::size_t>(c)] = 2;
      }
      std::vector<int> halo;
      halo.reserve(owned.size());
      auto paint = [&](int nc) {
        if (mark[static_cast<std::size_t>(nc)] != 0) {
          return;
        }
        mark[static_cast<std::size_t>(nc)] = 1;
        if (cellStart[nc] != cellStart[nc + 1]) {
          halo.push_back(nc);
        }
      };
      for (int c : owned) {
        const int cz = c % ncz;
        const int cy = (c / ncz) % ncy;
        const int cx = c / nyz;
        for (int ox = -rx; ox <= rx; ++ox) {
          int nx = cx + ox;
          while (nx < 0) {
            nx += ncx;
          }
          while (nx >= ncx) {
            nx -= ncx;
          }
          for (int oy = -ry; oy <= ry; ++oy) {
            int ny = cy + oy;
            while (ny < 0) {
              ny += ncy;
            }
            while (ny >= ncy) {
              ny -= ncy;
            }
            const int row = nx * ncy + ny;
            for (int oz = -rz; oz <= rz; ++oz) {
              if (ox == 0 && oy == 0 && oz == 0) {
                continue;
              }
              int nz = cz + oz;
              while (nz < 0) {
                nz += ncz;
              }
              while (nz >= ncz) {
                nz -= ncz;
              }
              paint(row * ncz + nz);
            }
          }
        }
      }
      std::vector<int> begin(static_cast<std::size_t>(nCells), 0);
      std::vector<int> end(static_cast<std::size_t>(nCells), 0);
      int total = 0;
      auto take = [&](int c) {
        const int n = cellStart[c + 1] - cellStart[c];
        begin[static_cast<std::size_t>(c)] = total;
        total += n;
        end[static_cast<std::size_t>(c)] = total;
      };
      for (int c : owned) {
        take(c);
      }
      for (int c : halo) {
        take(c);
      }
      std::vector<double> lx(static_cast<std::size_t>(total));
      std::vector<double> ly(static_cast<std::size_t>(total));
      std::vector<double> lz(static_cast<std::size_t>(total));
      std::vector<int> lid(static_cast<std::size_t>(total));
      int cursor = 0;
      auto copyCell = [&](int c) {
        for (int s = cellStart[c]; s < cellStart[c + 1]; ++s, ++cursor) {
          const std::size_t u = static_cast<std::size_t>(cursor);
          lx[u] = x[s];
          ly[u] = y[s];
          lz[u] = z[s];
          lid[u] = ids[s];
        }
      };
      for (int c : owned) {
        copyCell(c);
      }
      for (int c : halo) {
        copyCell(c);
      }
      const double *xs = lx.data();
      const double *ys = ly.data();
      const double *zs = lz.data();
      const int *id = lid.data();
      thread_local std::vector<PairDir> tlsBag;
      auto &bag = directRows == nullptr ? bags[static_cast<std::size_t>(tid)] : tlsBag;
      bag.clear();
      if (bag.capacity() < static_cast<std::size_t>(std::max(total, 1)) * 4) {
        bag.reserve(static_cast<std::size_t>(std::max(total, 1)) * 4);
      }
      int *mine = directRows == nullptr
                      ? deg.data() + static_cast<std::size_t>(tid) * static_cast<std::size_t>(nop)
                      : nullptr;
      const auto search0 = Clock::now();
      struct Nbr {
        int begin;
        int end;
        double sx, sy, sz;
        // 0 same cell (larger local slot writes both rows), 1 owned forward
        // cell (this side writes both rows), 2 ghost (this side writes its row).
        char rel;
      };
      std::vector<Nbr> wide;
      for (int c : owned) {
        const int b = begin[static_cast<std::size_t>(c)];
        const int e = end[static_cast<std::size_t>(c)];
        const int cz = c % ncz;
        const int cy = (c / ncz) % ncy;
        const int cx = c / nyz;
        Nbr nbr[27];
        int nn = 0;
        if (near) {
          for (int ox = -1; ox <= 1; ++ox) {
            int nx = cx + ox;
            int ix = 0;
            if (nx >= ncx) {
              nx -= ncx;
              ix = 1;
            } else if (nx < 0) {
              nx += ncx;
              ix = -1;
            }
            for (int oy = -1; oy <= 1; ++oy) {
              int ny = cy + oy;
              int iy = 0;
              if (ny >= ncy) {
                ny -= ncy;
                iy = 1;
              } else if (ny < 0) {
                ny += ncy;
                iy = -1;
              }
              const int row = nx * ncy + ny;
              for (int oz = -1; oz <= 1; ++oz) {
                int nz = cz + oz;
                int iz = 0;
                if (nz >= ncz) {
                  nz -= ncz;
                  iz = 1;
                } else if (nz < 0) {
                  nz += ncz;
                  iz = -1;
                }
                const int nc = row * ncz + nz;
                const int nb = begin[static_cast<std::size_t>(nc)];
                const int ne = end[static_cast<std::size_t>(nc)];
                if (nb == ne) {
                  continue;
                }
                const bool ownedNbr = mark[static_cast<std::size_t>(nc)] == 2;
                const bool backward = oz < 0 || (oz == 0 && oy < 0) || (oz == 0 && oy == 0 && ox < 0);
                if (ownedNbr && backward) {
                  continue;
                }
                const int k = (ix + 1) * 9 + (iy + 1) * 3 + (iz + 1);
                char rel = 2;
                if (ownedNbr) {
                  rel = (ox == 0 && oy == 0 && oz == 0) ? 0 : 1;
                }
                nbr[nn++] = {nb, ne, shx[k], shy[k], shz[k], rel};
              }
            }
          }
        } else {
          wide.clear();
          for (int ox = -rx; ox <= rx; ++ox) {
            int nx = cx + ox;
            int ix = 0;
            while (nx < 0) {
              nx += ncx;
              --ix;
            }
            while (nx >= ncx) {
              nx -= ncx;
              ++ix;
            }
            for (int oy = -ry; oy <= ry; ++oy) {
              int ny = cy + oy;
              int iy = 0;
              while (ny < 0) {
                ny += ncy;
                --iy;
              }
              while (ny >= ncy) {
                ny -= ncy;
                ++iy;
              }
              for (int oz = -rz; oz <= rz; ++oz) {
                int nz = cz + oz;
                int iz = 0;
                while (nz < 0) {
                  nz += ncz;
                  --iz;
                }
                while (nz >= ncz) {
                  nz -= ncz;
                  ++iz;
                }
                const int nc = (nx * ncy + ny) * ncz + nz;
                const int nb = begin[static_cast<std::size_t>(nc)];
                const int ne = end[static_cast<std::size_t>(nc)];
                if (nb == ne) {
                  continue;
                }
                wide.push_back({nb, ne, static_cast<double>(ix) * len[0],
                                static_cast<double>(iy) * len[1], static_cast<double>(iz) * len[2],
                                2});
                nn = static_cast<int>(wide.size());
              }
            }
          }
        }
        for (int i = b; i < e; ++i) {
          const int ii = id[static_cast<std::size_t>(i)];
          const double xi = xs[static_cast<std::size_t>(i)];
          const double yi = ys[static_cast<std::size_t>(i)];
          const double zi = zs[static_cast<std::size_t>(i)];
          const Nbr *span = near ? nbr : wide.data();
          const int nspan = near ? nn : static_cast<int>(wide.size());
          for (int t = 0; t < nspan; ++t) {
            const Nbr nb = span[t];
            for (int j = nb.begin; j < nb.end; ++j) {
              if (nb.rel == 0 && j <= i) {
                continue;
              }
              const double dxv = xi - xs[static_cast<std::size_t>(j)] - nb.sx;
              const double dyv = yi - ys[static_cast<std::size_t>(j)] - nb.sy;
              const double dzv = zi - zs[static_cast<std::size_t>(j)] - nb.sz;
              if (dxv * dxv + dyv * dyv + dzv * dzv >= rc2) {
                continue;
              }
              const int jj = id[static_cast<std::size_t>(j)];
              if (jj == ii) {
                continue;
              }
              if (nb.rel == 2) {
                bag.push_back({ii, jj});
                if (mine != nullptr) {
                  mine[ii]++;
                }
              } else {
                bag.push_back({ii, jj});
                bag.push_back({jj, ii});
                if (mine != nullptr) {
                  mine[ii]++;
                  mine[jj]++;
                }
              }
            }
          }
        }
      }
      if (directRows != nullptr) {
        const auto commit0 = Clock::now();
        std::vector<int> start(static_cast<std::size_t>(nop) + 1, 0);
        for (const PairDir &p : bag) {
          start[static_cast<std::size_t>(p.a) + 1]++;
        }
        for (int i = 0; i < nop; ++i) {
          start[static_cast<std::size_t>(i) + 1] += start[static_cast<std::size_t>(i)];
        }
        std::vector<int> neigh(bag.size());
        std::vector<int> cursor = start;
        for (const PairDir &p : bag) {
          neigh[static_cast<std::size_t>(cursor[static_cast<std::size_t>(p.a)]++)] = p.b;
        }
        for (int c : owned) {
          for (int s = begin[static_cast<std::size_t>(c)]; s < end[static_cast<std::size_t>(c)]; ++s) {
            const int ii = id[static_cast<std::size_t>(s)];
            auto &row = (*directRows)[static_cast<std::size_t>(ii)];
            if (indexToID != nullptr && indexToID[static_cast<std::size_t>(ii)] == -1) {
              row.clear();
              continue;
            }
            const int a = start[static_cast<std::size_t>(ii)];
            const int b = start[static_cast<std::size_t>(ii) + 1];
            if (b - a > 1) {
              sortRow(neigh.data() + a, b - a);
            }
            thread_local std::vector<int> tmp;
            tmp.clear();
            int prev = -1;
            for (int k = a; k < b; ++k) {
              const int jj = neigh[static_cast<std::size_t>(k)];
              if (jj == prev) {
                continue;
              }
              prev = jj;
              if (indexToID != nullptr && indexToID[static_cast<std::size_t>(jj)] == -1) {
                continue;
              }
              tmp.push_back(indexToID != nullptr ? indexToID[static_cast<std::size_t>(jj)] : jj);
            }
            const int nwrite = static_cast<int>(tmp.size()) + 1;
            if (static_cast<int>(row.capacity()) < nwrite) {
              row.reserve(static_cast<std::size_t>(nwrite));
            }
            row.resize(static_cast<std::size_t>(nwrite));
            row[0] = indexToID != nullptr ? indexToID[static_cast<std::size_t>(ii)] : ii;
            for (int k = 0; k < static_cast<int>(tmp.size()); ++k) {
              row[static_cast<std::size_t>(k + 1)] = tmp[static_cast<std::size_t>(k)];
            }
          }
        }
        if (profileOn() && tid == 0) {
          std::fprintf(stderr, "commit_ms %.3f\n", msSince(commit0));
        }
        bag.clear();
      }
      if (profileOn() && tid == 0) {
        std::fprintf(stderr,
                     "block_ms %.3f search_ms %.3f owned_cells %zu halo_cells %zu local %d\n",
                     msSince(block0) - msSince(search0), msSince(search0), owned.size(), halo.size(),
                     total);
      }
    }
  }
  return badTeam == 0;
}
#endif

bool orthoFullCsr(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<int> *subset, double rcutoff, std::vector<int> &offsets,
    std::vector<int> &cols, bool half, std::vector<std::vector<int>> *directRows,
    const int *indexToID, bool *wroteRows) {
  const int nop = yCloud.nop;
  if (nop < 0 || !(rcutoff > 0.0)) {
    return false;
  }
  if (yCloud.box.size() >= 6) {
    return tiltFullCsr(yCloud, subset, rcutoff, offsets, cols, half);
  }
  if (yCloud.box.size() != 3) {
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
  // Eight threads lose to one thread below this size. A 1000-molecule
  // ring census was slower with a team than without one.
  struct SerialBelow {
    int saved = 1;
    bool on = false;
    explicit SerialBelow(int n) {
#ifdef SEAMS_HAS_OPENMP
      if (n < 2048 && !omp_in_parallel()) {
        saved = std::max(1, omp_get_max_threads());
        if (saved > 1) {
          omp_set_num_threads(1);
          on = true;
        }
      }
#else
      (void)n;
#endif
    }
    ~SerialBelow() {
#ifdef SEAMS_HAS_OPENMP
      if (on) {
        omp_set_num_threads(saved);
      }
#endif
    }
  } serialBelow(nSubset);
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
  // Same bound as the tilted list: two cells keep a single image.
  const bool fast = reach[0] == 1 && reach[1] == 1 && reach[2] == 1 &&
                    ncell[0] >= 2 && ncell[1] >= 2 && ncell[2] >= 2;
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
#ifdef SEAMS_HAS_LINKCELL
  // An unsized Rayon pool keeps a worker on every cpu. Size it to this
  // OpenMP width before the stencil so those workers are not a second team.
  // Inside a live team the pool stays at one: the team already owns the cpus.
  lc_set_threads(nested ? 1 : maxThreads);
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
  double shx[27], shy[27], shz[27];
  fillOrthoShift(len[0], len[1], len[2], shx, shy, shz);

  const auto tCount = Clock::now();
  const bool stashOk = fast && !half;
  if (stashOk &&
      sc.stash.size() < static_cast<std::size_t>(nop) * static_cast<std::size_t>(kStashCap)) {
    sc.stash.resize(static_cast<std::size_t>(nop) * static_cast<std::size_t>(kStashCap));
  }
  // One half-stencil walk. Each unordered pair is tested once and written
  // into both rows. Private degree rows keep the threads off one counts line.
  //
  // The linkcell block copies each thread's cells and halo before the
  // stencil. On ice that copy is slower than the one-pass count walk, and
  // the block never sized the Rayon pool, so eight threads lost to one
  // thread and to vesin. The CSR and the row-vector callers both take the
  // count walk. The caller builds row vectors from that CSR.
  if (fast) {
    if (half) {
      passFast<false, true>(x, y, z, sc.ids.data(), sc.cellStart.data(), ncell[0], ncell[1],
                            ncell[2], shx, shy, shz, rc2, sc.counts.data(), nullptr, nullptr);
    } else if (maxThreads == 1) {
      passFastStashHalf(x, y, z, sc.ids.data(), sc.cellStart.data(), ncell[0], ncell[1],
                        ncell[2], shx, shy, shz, rc2, sc.counts.data(), sc.stash.data(),
                        &sc.overflowed);
    } else {
      passFastStash(x, y, z, sc.ids.data(), sc.cellStart.data(), ncell[0], ncell[1], ncell[2],
                    shx, shy, shz, rc2, sc.counts.data(), sc.stash.data(), &sc.overflowed);
    }
  } else {
    passSafe<false>(x, y, z, sc.ids.data(), sc.cellStart.data(), ncell[0], ncell[1],
                    ncell[2], reach[0], reach[1], reach[2], len[0], len[1], len[2], rc2,
                    sc.counts.data(), nullptr, nullptr, half);
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
    if (stashOk && sc.overflowed == 0) {
      scatterStash(sc.stash.data(), sc.counts.data(), offsets.data(), cols.data(), nop);
    } else if (fast) {
      if (half) {
        passFast<true, true>(x, y, z, sc.ids.data(), sc.cellStart.data(), ncell[0], ncell[1],
                             ncell[2], shx, shy, shz, rc2, nullptr, offsets.data(), cols.data());
      } else {
        passFast<true, false>(x, y, z, sc.ids.data(), sc.cellStart.data(), ncell[0], ncell[1],
                              ncell[2], shx, shy, shz, rc2, nullptr, offsets.data(), cols.data());
      }
    } else {
      passSafe<true>(x, y, z, sc.ids.data(), sc.cellStart.data(), ncell[0], ncell[1],
                     ncell[2], reach[0], reach[1], reach[2], len[0], len[1], len[2], rc2,
                     nullptr, offsets.data(), cols.data(), half);
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

nneigh::CutoffSkin::CutoffSkin(double cutoff, double skin)
    : cutoff_(cutoff), skin_(skin) {}

bool nneigh::CutoffSkin::moved(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud) const {
  if (!have_ || nop_ != yCloud.nop || static_cast<int>(x0_.size()) != yCloud.nop) {
    return true;
  }
  if (yCloud.box.size() != box0_.size() || yCloud.boxLow.size() != boxLow0_.size()) {
    return true;
  }
  for (std::size_t k = 0; k < yCloud.boxLow.size(); ++k) {
    if (std::fabs(yCloud.boxLow[k] - boxLow0_[k]) > 0.5 * skin_) {
      return true;
    }
  }
  const double trigger = 0.5 * skin_;
  for (std::size_t k = 0; k < yCloud.box.size(); ++k) {
    if (std::fabs(yCloud.box[k] - box0_[k]) > trigger) {
      return true;
    }
  }
  const gen::FracBox frame = gen::makeFracBox(yCloud);
  const double triggerSq = trigger * trigger;
  for (int i = 0; i < yCloud.nop; ++i) {
    if (subset_ != nullptr) {
      // A subset change is handled by the caller comparing pointers.
    }
    const auto &p = yCloud.pts[static_cast<std::size_t>(i)];
    const auto dr = gen::fracDelta(frame, p.x, p.y, p.z, x0_[static_cast<std::size_t>(i)],
                                   y0_[static_cast<std::size_t>(i)],
                                   z0_[static_cast<std::size_t>(i)]);
    if (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] > triggerSq) {
      return true;
    }
  }
  return false;
}

void nneigh::CutoffSkin::filter(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  const int nop = yCloud.nop;
  const double rc2 = cutoff_ * cutoff_;
  const gen::FracBox frame = gen::makeFracBox(yCloud);
  std::vector<int> counts(static_cast<std::size_t>(nop), 0);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nop >= 256)
#endif
  for (int i = 0; i < nop; ++i) {
    int n = 0;
    const int a = wideOff_[static_cast<std::size_t>(i)];
    const int b = wideOff_[static_cast<std::size_t>(i) + 1];
    for (int k = a; k < b; ++k) {
      const int j = wideCols_[static_cast<std::size_t>(k)];
      if (j < 0 || j >= nop || j == i) {
        continue;
      }
      if (half_ && j < i) {
        continue;
      }
      if (gen::periodicDistSq(frame, yCloud, i, j) < rc2) {
        ++n;
      }
    }
    counts[static_cast<std::size_t>(i)] = n;
  }
  offsets_.assign(static_cast<std::size_t>(nop) + 1, 0);
  for (int i = 0; i < nop; ++i) {
    offsets_[static_cast<std::size_t>(i) + 1] =
        offsets_[static_cast<std::size_t>(i)] + counts[static_cast<std::size_t>(i)];
  }
  cols_.resize(static_cast<std::size_t>(offsets_[static_cast<std::size_t>(nop)]));
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nop >= 256)
#endif
  for (int i = 0; i < nop; ++i) {
    int slot = offsets_[static_cast<std::size_t>(i)];
    const int a = wideOff_[static_cast<std::size_t>(i)];
    const int b = wideOff_[static_cast<std::size_t>(i) + 1];
    for (int k = a; k < b; ++k) {
      const int j = wideCols_[static_cast<std::size_t>(k)];
      if (j < 0 || j >= nop || j == i || (half_ && j < i)) {
        continue;
      }
      if (gen::periodicDistSq(frame, yCloud, i, j) < rc2) {
        cols_[static_cast<std::size_t>(slot++)] = j;
      }
    }
    const int n = slot - offsets_[static_cast<std::size_t>(i)];
    std::sort(cols_.begin() + offsets_[static_cast<std::size_t>(i)],
              cols_.begin() + offsets_[static_cast<std::size_t>(i)] + n);
  }
}

bool nneigh::CutoffSkin::update(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<int> *subset, bool half) {
  const bool sameSubset = subset_ == subset;
  const bool sameHalf = !have_ || half_ == half;
  half_ = half;
  reused_ = have_ && sameSubset && sameHalf && !moved(yCloud);
  if (!reused_) {
    if (!orthoFullCsr(yCloud, subset, cutoff_ + skin_, wideOff_, wideCols_, half_)) {
      have_ = false;
      offsets_.clear();
      cols_.clear();
      return false;
    }
    x0_.resize(static_cast<std::size_t>(yCloud.nop));
    y0_.resize(static_cast<std::size_t>(yCloud.nop));
    z0_.resize(static_cast<std::size_t>(yCloud.nop));
    for (int i = 0; i < yCloud.nop; ++i) {
      x0_[static_cast<std::size_t>(i)] = yCloud.pts[static_cast<std::size_t>(i)].x;
      y0_[static_cast<std::size_t>(i)] = yCloud.pts[static_cast<std::size_t>(i)].y;
      z0_[static_cast<std::size_t>(i)] = yCloud.pts[static_cast<std::size_t>(i)].z;
    }
    box0_ = yCloud.box;
    boxLow0_ = yCloud.boxLow;
    nop_ = yCloud.nop;
    subset_ = subset;
    have_ = true;
  }
  filter(yCloud);
  return true;
}

} // namespace nneigh
