// Host cageAffiliation versus the device (offload) cage pipeline, plus a
// brute-force reference of the order-invariant HC/DDC predicates.
//
// usage: check_device_host <lammpstrj> <frame> <type> <cutoff> [nRelabel]
//                          [cutoff|knn]
//        check_device_host levi:<q> [nRelabel]
//
// levi:<q> is the point-line incidence graph of the projective plane
// PG(2, q) for a prime q: (q+1)-regular and of girth 6, so every 6-cycle
// is a primitive ring, with q^3 (q+1) / 2 of them through each vertex
// (12, 54 and 375 for q = 2, 3 and 5).
//
// Prints one line per method with six-ring, HC/DDC ring and HC/DDC atom
// counts. The device path runs from its default buffers and from buffers
// of one entry, so every pass overflows until the buffers have grown.
// nRelabel > 0 reruns host and device on randomly relabelled copies of
// the graph; per atom, the labels may not change. Exits 1 on any
// disagreement in counts or per-atom labels.

#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>
#include <tum_offload.hpp>

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <random>
#include <string>
#include <vector>

namespace {

using Graph = std::vector<std::vector<int>>;

struct Labels {
  int nSix = 0;
  int hcRings = 0;
  int ddcRings = 0;
  int hcAtoms = 0;
  int ddcAtoms = 0;
  std::vector<int> atomHc;
  std::vector<int> atomDdc;
};

class Reference {
public:
  Reference(const Graph &rings, const Graph &nList)
      : rings_(rings), nList_(nList) {
    through_.resize(nList.size());
    for (int r = 0; r < static_cast<int>(rings.size()); ++r) {
      for (const int a : rings[static_cast<std::size_t>(r)]) {
        through_[static_cast<std::size_t>(a)].push_back(r);
      }
    }
  }

  Labels run() const {
    const int nR = static_cast<int>(rings_.size());
    std::vector<char> hc(static_cast<std::size_t>(nR), 0);
    int left = nR;
    for (int i = 0; i < nR && left > 0; ++i) {
      std::vector<int> cand;
      for (const int a : rings_[static_cast<std::size_t>(i)]) {
        const auto &row = nList_[static_cast<std::size_t>(a)];
        for (std::size_t n = 1; n < row.size(); ++n) {
          for (const int s : through_[static_cast<std::size_t>(row[n])]) {
            cand.push_back(s);
          }
        }
      }
      std::sort(cand.begin(), cand.end());
      cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
      auto mark = [&](int m) {
        if (m < 0 || m >= nR || hc[static_cast<std::size_t>(m)]) {
          return;
        }
        hc[static_cast<std::size_t>(m)] = 1;
        --left;
      };
      for (const int j : cand) {
        if (left == 0) {
          break;
        }
        if (j == i || !beta(rings_[static_cast<std::size_t>(i)],
                            rings_[static_cast<std::size_t>(j)])) {
          continue;
        }
        mark(i);
        mark(j);
        for (const int k : prismatic(i, j)) {
          mark(k);
        }
      }
    }
    std::vector<char> ddc(static_cast<std::size_t>(nR), 0);
    for (int i = 0; i < nR; ++i) {
      if (hc[static_cast<std::size_t>(i)]) {
        continue;
      }
      std::vector<int> peri;
      if (epsilon(i, peri)) {
        ddc[static_cast<std::size_t>(i)] = 1;
        for (const int p : peri) {
          ddc[static_cast<std::size_t>(p)] = 1;
        }
      }
    }
    Labels c;
    c.nSix = nR;
    c.atomHc.assign(nList_.size(), 0);
    c.atomDdc.assign(nList_.size(), 0);
    for (std::size_t r = 0; r < rings_.size(); ++r) {
      c.hcRings += hc[r];
      c.ddcRings += ddc[r];
      for (const int a : rings_[r]) {
        c.atomHc[static_cast<std::size_t>(a)] |= hc[r];
        c.atomDdc[static_cast<std::size_t>(a)] |= ddc[r];
      }
    }
    c.hcAtoms =
        static_cast<int>(std::count(c.atomHc.begin(), c.atomHc.end(), 1));
    c.ddcAtoms =
        static_cast<int>(std::count(c.atomDdc.begin(), c.atomDdc.end(), 1));
    return c;
  }

private:
  const Graph &rings_;
  const Graph &nList_;
  Graph through_;

  bool bonded(int a, int b) const {
    const auto &row = nList_[static_cast<std::size_t>(a)];
    return std::find(row.begin() + 1, row.end(), b) != row.end();
  }

  static bool in(const std::vector<int> &r, int a) {
    return std::find(r.begin(), r.end(), a) != r.end();
  }

  static int common(const std::vector<int> &a, const std::vector<int> &b) {
    int n = 0;
    for (const int x : a) {
      n += in(b, x) ? 1 : 0;
    }
    return n;
  }

  // beta(R, R'): all 2 x 2 x 6 alignments, nothing borrowed from the device
  bool beta(const std::vector<int> &R, const std::vector<int> &S) const {
    if (common(R, S) != 0) {
      return false;
    }
    const std::size_t perms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2},
                                     {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (std::size_t p = 0; p < 2; ++p) {
      for (std::size_t q = 0; q < 2; ++q) {
        bool oddFree = true;
        for (std::size_t t = 0; t < 3 && oddFree; ++t) {
          for (const int x : R) {
            if (bonded(x, S[(q + 1 + 2 * t) % 6])) {
              oddFree = false;
              break;
            }
          }
        }
        if (!oddFree) {
          continue;
        }
        for (const auto &s : perms) {
          bool ok = true;
          for (std::size_t t = 0; t < 3 && ok; ++t) {
            ok = bonded(R[(p + 2 * t) % 6], S[(q + 2 * s[t]) % 6]);
          }
          if (ok) {
            return true;
          }
        }
      }
    }
    return false;
  }

  std::vector<int> ringsWith(const std::vector<int> &atoms, int skipA,
                             int skipB) const {
    std::vector<int> out;
    for (const int r : through_[static_cast<std::size_t>(atoms[0])]) {
      if (r == skipA || r == skipB) {
        continue;
      }
      bool all = true;
      for (const int a : atoms) {
        all = all && in(rings_[static_cast<std::size_t>(r)], a);
      }
      if (all) {
        out.push_back(r);
      }
    }
    return out;
  }

  std::vector<int> prismatic(int i, int j) const {
    std::vector<int> out;
    const auto &R = rings_[static_cast<std::size_t>(i)];
    const auto &S = rings_[static_cast<std::size_t>(j)];
    for (std::size_t q = 0; q < 6; ++q) {
      const std::vector<int> trip{R[q], R[(q + 1) % 6], R[(q + 2) % 6]};
      for (const int k : ringsWith(trip, i, j)) {
        int rest = 0;
        int inS = 0;
        for (const int a : rings_[static_cast<std::size_t>(k)]) {
          if (!in(trip, a)) {
            ++rest;
            inS += in(S, a) ? 1 : 0;
          }
        }
        if (rest == 3 && inS == 3) {
          out.push_back(k);
        }
      }
    }
    return out;
  }

  bool epsilon(int i, std::vector<int> &peri) const {
    const auto &R = rings_[static_cast<std::size_t>(i)];
    for (const int a : R) {
      if (static_cast<int>(through_[static_cast<std::size_t>(a)].size()) - 1 <
          3) {
        return false;
      }
    }
    Graph T(6);
    for (std::size_t k = 0; k < 6; ++k) {
      T[k] = ringsWith({R[k], R[(k + 1) % 6], R[(k + 2) % 6]}, i, -1);
      if (T[k].empty()) {
        return false;
      }
    }
    auto shareOne = [&](int a, int b, int c) {
      for (const int x : rings_[static_cast<std::size_t>(a)]) {
        if (in(rings_[static_cast<std::size_t>(b)], x) &&
            in(rings_[static_cast<std::size_t>(c)], x)) {
          return true;
        }
      }
      return false;
    };
    // Every choice of P_0..P_5, the even and odd triples filtered as
    // soon as they are chosen
    bool any = false;
    std::array<int, 6> P{};
    for (const int p0 : T[0])
      for (const int p2 : T[2])
        for (const int p4 : T[4]) {
          if (p0 == p2 || p0 == p4 || p2 == p4 || !shareOne(p0, p2, p4)) {
            continue;
          }
          for (const int p1 : T[1])
            for (const int p3 : T[3])
              for (const int p5 : T[5]) {
                P = {p0, p1, p2, p3, p4, p5};
                std::array<int, 6> u = P;
                std::sort(u.begin(), u.end());
                if (std::unique(u.begin(), u.end()) != u.end()) {
                  continue;
                }
                if (!shareOne(p1, p3, p5)) {
                  continue;
                }
                // The pair conditions for some starting atom and direction
                // of the ring: triplet k of the relabelled ring is T[s(k)]
                bool aligned = false;
                for (int r = 0; r < 6 && !aligned; ++r) {
                  for (int flip = 0; flip < 2 && !aligned; ++flip) {
                    auto at = [&](int k) {
                      const int slot = flip ? (r - k + 12) % 6 : (r + k) % 6;
                      return static_cast<std::size_t>(
                          P[static_cast<std::size_t>(slot)]);
                    };
                    aligned = common(rings_[at(0)], rings_[at(2)]) >= 3 &&
                              common(rings_[at(1)], rings_[at(3)]) >= 3 &&
                              common(rings_[at(2)], rings_[at(4)]) >= 3 &&
                              common(rings_[at(3)], rings_[at(5)]) >= 3;
                  }
                }
                if (!aligned) {
                  continue;
                }
                any = true;
                peri.insert(peri.end(), P.begin(), P.end());
              }
        }
    std::sort(peri.begin(), peri.end());
    peri.erase(std::unique(peri.begin(), peri.end()), peri.end());
    return any;
  }
};

Labels fromCageCounts(const tum::CageCounts &c) {
  return {c.nSix,      c.nHcRings, c.nDdcRings, c.nHcAtoms,
          c.nDdcAtoms, c.atomHc,   c.atomDdc};
}

void print(const char *what, const Labels &c) {
  std::printf("%-22s nSix %6d  HC rings %6d  DDC rings %6d  HC atoms %6d  "
              "DDC atoms %6d\n",
              what, c.nSix, c.hcRings, c.ddcRings, c.hcAtoms, c.ddcAtoms);
}

void printPath(const char *what, const tum::CageCounts &c) {
  print(what, fromCageCounts(c));
  std::printf("%-22s retries %d  host fallback %s\n", "", c.retries,
              c.hostFallback ? "yes" : "no");
}

bool same(const Labels &a, const Labels &b) {
  return a.nSix == b.nSix && a.hcRings == b.hcRings &&
         a.ddcRings == b.ddcRings && a.hcAtoms == b.hcAtoms &&
         a.ddcAtoms == b.ddcAtoms && a.atomHc == b.atomHc &&
         a.atomDdc == b.atomDdc;
}

// Per-atom labels of a relabelled graph read back in the original labels
Labels unpermute(Labels c, const std::vector<int> &perm) {
  std::vector<int> hc(perm.size()), ddc(perm.size());
  for (std::size_t a = 0; a < perm.size(); ++a) {
    hc[a] = c.atomHc[static_cast<std::size_t>(perm[a])];
    ddc[a] = c.atomDdc[static_cast<std::size_t>(perm[a])];
  }
  c.atomHc = std::move(hc);
  c.atomDdc = std::move(ddc);
  return c;
}

Graph relabel(const Graph &nList, const std::vector<int> &perm) {
  // perm[old] = new; rows keep row[0] = self, neighbours shuffled too
  Graph out(nList.size());
  std::mt19937 rng(static_cast<unsigned>(perm[0]) * 7919U + 17U);
  for (std::size_t a = 0; a < nList.size(); ++a) {
    std::vector<int> row;
    row.push_back(perm[a]);
    for (std::size_t n = 1; n < nList[a].size(); ++n) {
      row.push_back(perm[static_cast<std::size_t>(nList[a][n])]);
    }
    std::shuffle(row.begin() + 1, row.end(), rng);
    out[static_cast<std::size_t>(perm[a])] = row;
  }
  return out;
}

// Points then lines of PG(2, q), q prime, as normalised vectors of GF(q)^3
Graph leviGraph(int q) {
  std::vector<std::array<int, 3>> v;
  for (int a = 0; a < q; ++a) {
    for (int b = 0; b < q; ++b) {
      for (int c = 0; c < q; ++c) {
        const int lead = a != 0 ? a : (b != 0 ? b : c);
        if (lead == 1) {
          v.push_back({a, b, c});
        }
      }
    }
  }
  const int m = static_cast<int>(v.size());
  Graph nList(static_cast<std::size_t>(2 * m));
  for (int i = 0; i < 2 * m; ++i) {
    nList[static_cast<std::size_t>(i)].push_back(i);
  }
  for (int p = 0; p < m; ++p) {
    const auto &vp = v[static_cast<std::size_t>(p)];
    for (int l = 0; l < m; ++l) {
      const auto &vl = v[static_cast<std::size_t>(l)];
      const int dot = vp[0] * vl[0] + vp[1] * vl[1] + vp[2] * vl[2];
      if (dot % q == 0) {
        nList[static_cast<std::size_t>(p)].push_back(m + l);
        nList[static_cast<std::size_t>(m + l)].push_back(p);
      }
    }
  }
  return nList;
}

bool isPrime(int q) {
  if (q < 2) {
    return false;
  }
  for (int d = 2; d * d <= q; ++d) {
    if (q % d == 0) {
      return false;
    }
  }
  return true;
}

} // namespace

int main(int argc, char **argv) {
  const std::string source = argc > 1 ? argv[1] : "";
  const bool levi = source.rfind("levi:", 0) == 0;
  if (!levi && argc < 5) {
    std::fprintf(stderr,
                 "usage: %s <lammpstrj> <frame> <type> <cutoff> [nRelabel] "
                 "[cutoff|knn]\n       %s levi:<q> [nRelabel]\n",
                 argv[0], argv[0]);
    return 2;
  }
  Graph nList;
  int nRelabel = 0;
  if (levi) {
    const int q = std::atoi(source.c_str() + 5);
    if (!isPrime(q)) {
      std::fprintf(stderr, "levi:<q> needs a prime q\n");
      return 2;
    }
    nRelabel = argc > 2 ? std::atoi(argv[2]) : 0;
    nList = leviGraph(q);
    std::printf("Levi graph of PG(2,%d) atoms %zu degree %d\n", q, nList.size(),
                q + 1);
  } else {
    const int frame = std::atoi(argv[2]);
    const int type = std::atoi(argv[3]);
    const double cutoff = std::atof(argv[4]);
    nRelabel = argc > 5 ? std::atoi(argv[5]) : 0;
    const std::string graph = argc > 6 ? argv[6] : "cutoff";
    molSys::PointCloud<molSys::Point<double>, double> cloud;
    cloud = sinp::readLammpsTrjO(source, frame, cloud, type);
    // knn: the mutual 4-NN graph the CUDA batch path builds
    auto nl =
        graph == "knn"
            ? nneigh::kNearestNeighbourList(cloud, 4, cutoff + 1.5, type, true)
            : nneigh::neighListO(cutoff, cloud, type);
    nList = nneigh::neighbourListByIndex(cloud, nl);
    std::printf("%s frame %d type %d cutoff %.2f graph %s atoms %d\n",
                source.c_str(), frame, type, cutoff, graph.c_str(), cloud.nop);
  }

  tum::SpecializedCaps one;
  one.rowCap = 1;
  one.ringCap = 1;
  one.pairCap = 1;
  tum::SpecializedCaps target;
  target.forceTarget = true;

  const auto host = tum::hostCageCounts(nList);
  const auto device = tum::specializedCageCounts(nList);
  const auto grown = tum::specializedCageCounts(nList, one);
  const auto forced = tum::specializedCageCounts(nList, target);
  Graph six;
  for (auto &r : primitive::ringNetwork(nList, 6)) {
    if (r.size() == 6) {
      six.push_back(r);
    }
  }
  const Labels ref = Reference(six, nList).run();
  {
    std::vector<int> per(nList.size(), 0);
    for (const auto &r : six) {
      for (const int a : r) {
        ++per[static_cast<std::size_t>(a)];
      }
    }
    const int mx = per.empty() ? 0 : *std::max_element(per.begin(), per.end());
    auto above = [&](int cap) {
      return static_cast<int>(std::count_if(per.begin(), per.end(),
                                            [cap](int v) { return v > cap; }));
    };
    std::printf("max six-rings per atom %d, atoms above 16: %d, above 32: %d\n",
                mx, above(16), above(32));
  }

  const Labels h = fromCageCounts(host);
  const Labels d = fromCageCounts(device);
  const Labels g = fromCageCounts(grown);
  const Labels f = fromCageCounts(forced);
  print("host cageAffiliation", h);
  printPath("device (offload)", device);
  printPath("device (caps 1)", grown);
  printPath("device (target)", forced);
  print("spec reference", ref);
  bool ok =
      same(d, ref) && same(h, ref) && same(h, d) && same(g, d) && same(f, d);
  std::printf("device==reference %s  host==reference %s  host==device %s  "
              "caps 1==device %s  target==device %s  (counts and per-atom)\n",
              same(d, ref) ? "yes" : "NO", same(h, ref) ? "yes" : "NO",
              same(h, d) ? "yes" : "NO", same(g, d) ? "yes" : "NO",
              same(f, d) ? "yes" : "NO");

  std::mt19937 rng(12345U);
  for (int t = 0; t < nRelabel; ++t) {
    std::vector<int> perm(nList.size());
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), rng);
    const auto pg = relabel(nList, perm);
    const Labels rd =
        unpermute(fromCageCounts(tum::specializedCageCounts(pg)), perm);
    const Labels rh = unpermute(fromCageCounts(tum::hostCageCounts(pg)), perm);
    char label[32];
    std::snprintf(label, sizeof(label), "relabel %d", t);
    print(label, rd);
    const bool agree = same(rd, d) && same(rh, h);
    std::printf("%-22s device==unrelabelled %s  host==unrelabelled %s\n", "",
                same(rd, d) ? "yes" : "NO", same(rh, h) ? "yes" : "NO");
    ok = ok && agree;
  }
  std::printf("%s\n", ok ? "PASS" : "FAIL");
  return ok ? 0 : 1;
}
