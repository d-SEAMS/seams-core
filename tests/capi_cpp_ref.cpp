/* C++ path for the C API suite. Reads one frame file and prints the
 * same keys the C and Python drivers print. The calls are the library
 * functions, not the C wrappers.
 */
#include <bond.hpp>
#include <bop.hpp>
#include <brine_front.hpp>
#include <cage.hpp>
#include <cage_affiliation.hpp>
#include <cage_enum.hpp>
#include <chill_plus_scheme.hpp>
#include <cluster.hpp>
#include <density.hpp>
#include <edge_set.hpp>
#include <franzblau.hpp>
#include <local_symmetry.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <nucleation_pathway.hpp>
#include <rdf.hpp>
#include <ring.hpp>
#include <site.hpp>
#include <topo_bulk.hpp>
#include <topo_fingerprint.hpp>
#include <topo_one_dim.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

struct Frame {
  int n = 0;
  std::vector<double> xyz;
  std::vector<int> type;
  double box[9] = {};
};

Frame load(const char *path) {
  std::ifstream in(path);
  if (!in) {
    std::fprintf(stderr, "cannot read %s\n", path);
    std::exit(1);
  }
  Frame frame;
  in >> frame.n;
  for (int i = 0; i < 9; ++i) {
    in >> frame.box[i];
  }
  frame.xyz.resize(static_cast<std::size_t>(frame.n) * 3);
  frame.type.resize(static_cast<std::size_t>(frame.n));
  for (int i = 0; i < frame.n; ++i) {
    in >> frame.xyz[static_cast<std::size_t>(3 * i)];
    in >> frame.xyz[static_cast<std::size_t>(3 * i + 1)];
    in >> frame.xyz[static_cast<std::size_t>(3 * i + 2)];
    in >> frame.type[static_cast<std::size_t>(i)];
  }
  if (!in) {
    std::fprintf(stderr, "short frame %s\n", path);
    std::exit(1);
  }
  return frame;
}

Cloud cloudOf(const Frame &frame) {
  Cloud cloud;
  cloud.nop = frame.n;
  cloud.currentFrame = 1;
  cloud.pts.resize(static_cast<std::size_t>(frame.n));
  for (int i = 0; i < frame.n; ++i) {
    auto &point = cloud.pts[static_cast<std::size_t>(i)];
    point.x = frame.xyz[static_cast<std::size_t>(3 * i)];
    point.y = frame.xyz[static_cast<std::size_t>(3 * i + 1)];
    point.z = frame.xyz[static_cast<std::size_t>(3 * i + 2)];
    point.type = frame.type[static_cast<std::size_t>(i)];
    point.atomID = i + 1;
    point.molID = point.atomID;
    cloud.idIndexMap[point.atomID] = i;
  }
  if (!nneigh::setBoxFromLammpsEdges(cloud.box, cloud.boxLow, frame.box[0], frame.box[1],
                                    frame.box[2], frame.box[3], frame.box[4], frame.box[5],
                                    frame.box[6], frame.box[7], frame.box[8])) {
    std::fprintf(stderr, "bad box\n");
    std::exit(1);
  }
  return cloud;
}

void ki(const char *key, int value) { std::printf("%s %d\n", key, value); }

void kd(const char *key, double value) { std::printf("%s %.17g\n", key, value); }

double fingerprint(const std::vector<double> &values) {
  long double sum = 0;
  for (int i = 0; i < static_cast<int>(values.size()); ++i) {
    double x = values[static_cast<std::size_t>(i)];
    if (!std::isfinite(x)) {
      x = 0;
    }
    sum += static_cast<long double>(x) * (i + 1);
  }
  return static_cast<double>(sum);
}

double meanOf(const std::vector<double> &values) {
  long double sum = 0;
  int n = 0;
  for (double x : values) {
    if (!std::isfinite(x)) {
      continue;
    }
    sum += x;
    ++n;
  }
  return n > 0 ? static_cast<double>(sum / n) : 0;
}

double medianOf(std::vector<double> values) {
  values.erase(std::remove_if(values.begin(), values.end(),
                              [](double x) { return !std::isfinite(x); }),
               values.end());
  if (values.empty()) {
    return 0;
  }
  const std::size_t mid = values.size() / 2;
  std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(mid),
                   values.end());
  if (values.size() % 2 == 1) {
    return values[mid];
  }
  const double upper = values[mid];
  std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(mid - 1),
                   values.end());
  return 0.5 * (values[mid - 1] + upper);
}

void hist9(const char *prefix, const Cloud &cloud) {
  int counts[9] = {};
  for (const auto &point : cloud.pts) {
    int label = static_cast<int>(point.iceType);
    if (label < 0 || label > 8) {
      label = 6;
    }
    counts[label] += 1;
  }
  for (int i = 0; i < 9; ++i) {
    std::printf("%s_%d %d\n", prefix, i, counts[i]);
  }
}

std::vector<std::vector<int>> indexRows(const Cloud &cloud,
                                        const std::vector<std::vector<int>> &rows) {
  auto idx = nneigh::neighbourListByIndex(cloud, rows);
  if (static_cast<int>(idx.size()) < cloud.nop) {
    idx.resize(static_cast<std::size_t>(cloud.nop));
  }
  for (int i = 0; i < cloud.nop; ++i) {
    auto &row = idx[static_cast<std::size_t>(i)];
    row.erase(std::remove(row.begin(), row.end(), i), row.end());
  }
  return idx;
}

bool crystal(int label) { return label == 0 || label == 1 || label == 7 || label == 8; }

Cloud unified(const Cloud &cloud, int type, int &typ) {
  Cloud work = cloud;
  if (type > 0) {
    typ = type;
    return work;
  }
  if (work.nop <= 0) {
    typ = 1;
    return work;
  }
  const int first = work.pts[0].type;
  bool mixed = first <= 0;
  for (int i = 1; i < work.nop && !mixed; ++i) {
    if (work.pts[static_cast<std::size_t>(i)].type != first) {
      mixed = true;
    }
  }
  if (mixed) {
    for (auto &point : work.pts) {
      point.type = 1;
    }
    typ = 1;
    return work;
  }
  typ = first > 0 ? first : 1;
  if (typ != first) {
    for (auto &point : work.pts) {
      point.type = typ;
    }
  }
  return work;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 3) {
    std::fprintf(stderr, "usage: capi_cpp_ref FRAME KIND\n");
    return 2;
  }
  const std::string kind = argv[2];
  const Frame frame = load(argv[1]);
  const double cutoff = 3.5;
  Cloud cloud = cloudOf(frame);
  ki("nop", cloud.nop);

  {
    auto idx = indexRows(cloud, nneigh::neighListO(cutoff, cloud, 1));
    int nnz = 0;
    int lo = cloud.nop == 0 ? 0 : 1000;
    int hi = 0;
    for (const auto &row : idx) {
      const int degree = static_cast<int>(row.size());
      nnz += degree;
      lo = std::min(lo, degree);
      hi = std::max(hi, degree);
    }
    if (cloud.nop == 0) {
      lo = 0;
    }
    ki("deg_min", lo);
    ki("deg_max", hi);
    ki("nnz", nnz);
  }
  {
    std::vector<int> offsets;
    std::vector<int> cols;
    if (!nneigh::orthoFullCsr(cloud, nullptr, cutoff, offsets, cols, false)) {
      std::fprintf(stderr, "csr failed\n");
      return 1;
    }
    ki("csr_nnz", static_cast<int>(cols.size()));
    const auto rings = primitive::ringNetwork(offsets, cols, 6);
    int hist[7] = {};
    for (const auto &ring : rings) {
      const int size = static_cast<int>(ring.size());
      if (size >= 0 && size < 7) {
        hist[size] += 1;
      }
    }
    for (int i = 0; i < 7; ++i) {
      std::printf("ring_%d %d\n", i, hist[i]);
    }
  }
  for (int order : {3, 4, 6, 8, 12}) {
    const auto result = chill::steinhardtQlCutoff(cloud, cutoff, order, 1, false);
    char key[32];
    std::snprintf(key, sizeof key, "q%d_mean", order);
    kd(key, meanOf(result.ql));
    std::snprintf(key, sizeof key, "q%d_fp", order);
    kd(key, fingerprint(result.ql));
  }
  {
    const auto result = chill::steinhardtQlCutoff(cloud, cutoff, 6, 1, true);
    kd("q6bar_mean", meanOf(result.qlBar));
    kd("q6bar_p50", medianOf(result.qlBar));
    kd("q6bar_fp", fingerprint(result.qlBar));
    int ice = 0;
    for (double q : result.qlBar) {
      if (std::isfinite(q) && q > 0.36) {
        ++ice;
      }
    }
    ki("q6bar_ice", ice);
  }
  {
    Cloud work = cloud;
    auto rows = nneigh::neighListO(cutoff, work, 1);
    chill::getCorrel(work, rows, false);
    chill::getIceTypeNoPrint(work, rows, false);
    hist9("chill", work);
  }
  {
    Cloud work = cloud;
    if (!chill::chillPlusPublishedCutoff(work, cutoff, 1)) {
      std::vector<int> offsets;
      std::vector<int> cols;
      if (!nneigh::orthoFullCsr(work, nullptr, cutoff, offsets, cols, false)) {
        std::fprintf(stderr, "published csr failed\n");
        return 1;
      }
      chill::chillPlusFromCsr(work, offsets, cols);
    }
    hist9("plus", work);
  }
  {
    Cloud work = cloud;
    auto rows = nneigh::kNearestNeighbourList(work, 4, cutoff + 1.5, 1, true);
    chill::applyChillPlus(work, rows, chill::ChillPlusRule::Ovito, false);
    hist9("ovito", work);
  }
  {
    const int n = cloud.nop;
    auto graphs = nneigh::kNearestNeighbourPair(cloud, 4, cutoff + 1.5, 1);
    std::vector<char> hc(static_cast<std::size_t>(n), 0);
    std::vector<char> ddc(static_cast<std::size_t>(n), 0);
    auto idxStrict = nneigh::neighbourListByIndex(cloud, graphs.first);
    auto idxUnion = nneigh::neighbourListByIndex(cloud, graphs.second);
    const bool same = graphs.first == graphs.second ||
                      seams::provenEqual(seams::edgeSetOf(idxStrict), seams::edgeSetOf(idxUnion));
    if (same) {
      auto six = ring::contractibleRings(primitive::ringNetwork(idxStrict, 6), cloud);
      std::vector<std::vector<int>> only6;
      for (const auto &ring : six) {
        if (ring.size() == 6) {
          only6.push_back(ring);
        }
      }
      const auto once = ring::cageAffiliation(only6, idxStrict);
      const auto aff = ring::seededCageAffiliation(only6, once, only6, once, idxStrict, false);
      const int m = std::min(n, static_cast<int>(aff.hc.size()));
      for (int i = 0; i < m; ++i) {
        hc[static_cast<std::size_t>(i)] = aff.hc[static_cast<std::size_t>(i)];
        ddc[static_cast<std::size_t>(i)] = aff.ddc[static_cast<std::size_t>(i)];
      }
    } else {
      auto sixStrict = ring::contractibleRings(primitive::ringNetwork(idxStrict, 6), cloud);
      auto sixUnion = ring::contractibleRings(primitive::ringNetwork(idxUnion, 6), cloud);
      std::vector<std::vector<int>> a;
      std::vector<std::vector<int>> b;
      for (const auto &ring : sixStrict) {
        if (ring.size() == 6) {
          a.push_back(ring);
        }
      }
      for (const auto &ring : sixUnion) {
        if (ring.size() == 6) {
          b.push_back(ring);
        }
      }
      const auto aff = ring::seededCageAffiliation(a, idxStrict, b, idxUnion, false);
      const int m = std::min(n, static_cast<int>(aff.hc.size()));
      for (int i = 0; i < m; ++i) {
        hc[static_cast<std::size_t>(i)] = aff.hc[static_cast<std::size_t>(i)];
        ddc[static_cast<std::size_t>(i)] = aff.ddc[static_cast<std::size_t>(i)];
      }
    }
    int nHex = 0;
    int nCubic = 0;
    int nWater = 0;
    int nBoth = 0;
    for (int i = 0; i < n; ++i) {
      const bool h = hc[static_cast<std::size_t>(i)] != 0;
      const bool d = ddc[static_cast<std::size_t>(i)] != 0;
      if (h && d) {
        ++nBoth;
      } else if (h) {
        ++nHex;
      } else if (d) {
        ++nCubic;
      } else {
        ++nWater;
      }
    }
    ki("cage_hex", nHex);
    ki("cage_cubic", nCubic);
    ki("cage_water", nWater);
    ki("cage_both", nBoth);
  }
  {
    std::vector<double> x(static_cast<std::size_t>(cloud.nop));
    std::vector<double> y(x.size());
    std::vector<double> z(x.size());
    for (int i = 0; i < cloud.nop; ++i) {
      x[static_cast<std::size_t>(i)] = cloud.pts[static_cast<std::size_t>(i)].x;
      y[static_cast<std::size_t>(i)] = cloud.pts[static_cast<std::size_t>(i)].y;
      z[static_cast<std::size_t>(i)] = cloud.pts[static_cast<std::size_t>(i)].z;
    }
    double lat[3][3] = {{frame.box[3], 0, 0}, {frame.box[6], frame.box[4], 0},
                        {frame.box[7], frame.box[8], frame.box[5]}};
    std::vector<int> iceKind(static_cast<std::size_t>(cloud.nop), -1);
    symm::classifyIceArrays(x.data(), y.data(), z.data(), cloud.nop, lat, iceKind.data(), 0.12);
    int ih = 0;
    int ic = 0;
    int isd = 0;
    int other = 0;
    for (int label : iceKind) {
      if (label == 0) {
        ++ih;
      } else if (label == 1) {
        ++ic;
      } else if (label == 2) {
        ++isd;
      } else {
        ++other;
      }
    }
    ki("sym_ih", ih);
    ki("sym_ic", ic);
    ki("sym_isd", isd);
    ki("sym_other", other);
  }
  {
    auto ids = nneigh::neighListLinkcell(cutoff, cloud, 1);
    if (ids.empty()) {
      ids = nneigh::neighListO(cutoff, cloud, 1);
    }
    auto idx = nneigh::neighbourListByIndex(cloud, ids);
    const auto signature = cage::Signature::parse("512");
    const int depth = std::max(signature.maxRingSize(), 3);
    const auto rings = ring::cageRings(primitive::ringNetwork(idx, depth), cloud, cutoff);
    const auto window = cage::parseEdgeWindow("");
    std::vector<std::array<double, 3>> xyz;
    const std::vector<std::array<double, 3>> *xyzPtr = nullptr;
    const std::vector<double> *boxPtr = nullptr;
    if (window.enabled) {
      xyz.resize(static_cast<std::size_t>(cloud.nop));
      for (int i = 0; i < cloud.nop; ++i) {
        const auto &point = cloud.pts[static_cast<std::size_t>(i)];
        xyz[static_cast<std::size_t>(i)] = {point.x, point.y, point.z};
      }
      xyzPtr = &xyz;
      boxPtr = &cloud.box;
    }
    const auto foundOf = cage::findBySignatures(rings, idx, std::vector<cage::Signature>{signature},
                                               {}, xyzPtr, boxPtr, window);
    ki("sig_cages", static_cast<int>(foundOf.at(0).size()));
    std::vector<char> seen(static_cast<std::size_t>(cloud.nop), 0);
    int atoms = 0;
    for (const auto &cageUnit : foundOf.at(0)) {
      for (int vertex : cageUnit.vertices) {
        if (vertex >= 0 && vertex < cloud.nop &&
            seen[static_cast<std::size_t>(vertex)] == 0) {
          seen[static_cast<std::size_t>(vertex)] = 1;
          ++atoms;
        }
      }
    }
    ki("sig_atoms", atoms);
  }
  {
    int typ = 1;
    Cloud work = unified(cloud, 0, typ);
    auto rows = nneigh::kNearestNeighbourList(work, 4, cutoff + 1.5, typ, true);
    auto idx = nneigh::neighbourListByIndex(work, rows);
    int nnz = 0;
    const int n = std::min(work.nop, static_cast<int>(idx.size()));
    for (int i = 0; i < n; ++i) {
      for (int id : idx[static_cast<std::size_t>(i)]) {
        if (id != i) {
          ++nnz;
        }
      }
    }
    ki("knn_nnz", nnz);
    const auto gr = rdf::partialRdf(work, typ, typ, 10.0, 20);
    kd("rdf_vol", gr.volume);
    long double gsum = 0;
    int pairs = 0;
    for (std::size_t i = 0; i < gr.g.size() && i < gr.count.size(); ++i) {
      gsum += gr.g[i];
      pairs += gr.count[i];
    }
    kd("rdf_gsum", static_cast<double>(gsum));
    ki("rdf_count", pairs);
    const double rho = gr.volume > 0.0 ? static_cast<double>(gr.nJ) / gr.volume : 0.0;
    const auto series = rdf::runningCN(gr, rho);
    kd("cn_last", series.empty() ? 0.0 : series.back());
    kd("cn", rdf::coordinationNumber(gr, cutoff, rho));
    const auto density = site::densityZ(cloud, 0, 10, 2);
    long double dens = 0;
    for (double value : density.rho) {
      dens += value;
    }
    kd("dens_sum", static_cast<double>(dens));
    Cloud heavy;
    Cloud hydro;
    heavy.box = cloud.box;
    heavy.boxLow = cloud.boxLow;
    hydro.box = cloud.box;
    hydro.boxLow = cloud.boxLow;
    for (int i = 0; i < cloud.nop; ++i) {
      const auto &point = cloud.pts[static_cast<std::size_t>(i)];
      if (point.type == 2) {
        hydro.pts.push_back(point);
      } else if (point.type == 1) {
        heavy.pts.push_back(point);
      }
    }
    heavy.nop = static_cast<int>(heavy.pts.size());
    hydro.nop = static_cast<int>(hydro.pts.size());
    for (int i = 0; i < heavy.nop; ++i) {
      heavy.idIndexMap[heavy.pts[static_cast<std::size_t>(i)].atomID] = i;
    }
    int bonds = 0;
    if (heavy.nop > 0 && hydro.nop > 0) {
      const int heavyType = heavy.pts[0].type;
      auto bondRows = nneigh::neighListO(cutoff, heavy, heavyType);
      auto net = bond::populateHbondsWithInputClouds(heavy, hydro, bondRows, 2.5, 30.0);
      int edges = 0;
      for (const auto &row : net) {
        if (row.size() > 1) {
          edges += static_cast<int>(row.size()) - 1;
        }
      }
      bonds = edges / 2;
    }
    ki("hbond_n", bonds);
    auto mutual = nneigh::kNearestNeighbourList(work, 4, cutoff + 1.5, typ, true);
    auto prismIdx = nneigh::neighbourListByIndex(work, mutual);
    const auto prismRings = primitive::ringNetwork(prismIdx, 6);
    int perfect = 0;
    int imperfect = 0;
    for (int size = 3; size <= 6; ++size) {
      std::vector<std::vector<int>> one;
      for (const auto &ring : prismRings) {
        if (static_cast<int>(ring.size()) == size) {
          one.push_back(ring);
        }
      }
      if (one.empty()) {
        continue;
      }
      std::vector<ring::strucType> kindOf(one.size());
      std::vector<double> rmsd(static_cast<std::size_t>(std::max(work.nop, 0)), -1.0);
      int nPerfect = 0;
      int nImperfect = 0;
      ring::findPrisms(one, kindOf, nPerfect, nImperfect, prismIdx, work, rmsd, false);
      perfect += nPerfect;
      imperfect += nImperfect;
    }
    ki("prism_p", perfect);
    ki("prism_i", imperfect);
    auto cutoffRows = nneigh::neighListO(cutoff, work, typ);
    auto monoIdx = nneigh::neighbourListByIndex(work, cutoffRows);
    const auto monoRings = primitive::ringNetwork(monoIdx, 6);
    int mono = 0;
    for (const auto &ring : monoRings) {
      if (ring.size() < 7) {
        ++mono;
      }
    }
    ki("mono_n", mono);
    const auto table = site::parseSiteSpec("1=waterO");
    int water = 0;
    for (const auto &point : cloud.pts) {
      if (table.of(point) == site::Kind::waterO) {
        ++water;
      }
    }
    ki("site_water", water);
    const auto ions = site::ionCloud(cloud, site::parseSiteSpec("1=cationHead,2=anion"));
    const auto contacts = nneigh::mutualNearestUnlike(ions, 1, 2);
    int cations = 0;
    int anions = 0;
    for (const auto &point : ions.pts) {
      if (point.type == 1) {
        ++cations;
      } else if (point.type == 2) {
        ++anions;
      }
    }
    ki("contact_n", static_cast<int>(contacts.size()));
    ki("contact_cat", cations);
    ki("contact_an", anions);
    const auto domainTable = site::parseSiteSpec("1=apolar");
    std::vector<bool> mask(static_cast<std::size_t>(cloud.nop), false);
    for (const int index : site::indicesOf(cloud, domainTable, site::Kind::apolar)) {
      if (index >= 0 && index < cloud.nop) {
        mask[static_cast<std::size_t>(index)] = true;
      }
    }
    const auto domainIdx = nneigh::getNewNeighbourListByIndex(cloud, cutoff);
    std::vector<std::vector<int>> ids(domainIdx.size());
    for (std::size_t i = 0; i < domainIdx.size(); ++i) {
      for (const int index : domainIdx[i]) {
        if (index >= 0 && index < cloud.nop) {
          ids[i].push_back(cloud.pts[static_cast<std::size_t>(index)].atomID);
        }
      }
    }
    const auto domain = clump::largestDomain(cloud, ids, mask);
    ki("domain_n", domain.subset);
    ki("domain_largest", domain.largest);
    kd("domain_p", domain.percolation);
    std::vector<int> colours;
    colours.reserve(static_cast<std::size_t>(work.nop));
    for (int i = 0; i < work.nop; ++i) {
      colours.push_back(work.pts[static_cast<std::size_t>(i)].type);
    }
    const auto fp = topo::fingerprint(monoIdx, 2, 7, colours);
    ki("fp_len", static_cast<int>(fp.key.size()));
    int census = 0;
    for (int count : fp.ringCensus) {
      census += count;
    }
    ki("fp_rings", census);
  }
  if (kind == "brine") {
    Cloud work = cloud;
    std::vector<int> water;
    std::vector<int> ions;
    for (int i = 0; i < work.nop; ++i) {
      if (work.pts[static_cast<std::size_t>(i)].type == 1) {
        water.push_back(i);
      } else if (work.pts[static_cast<std::size_t>(i)].type == 3) {
        ions.push_back(i);
      }
    }
    std::vector<int> labels(static_cast<std::size_t>(work.nop), 6);
    if (!water.empty()) {
      Cloud oxygen;
      oxygen.nop = static_cast<int>(water.size());
      oxygen.currentFrame = 1;
      oxygen.box = work.box;
      oxygen.boxLow = work.boxLow;
      oxygen.pts.resize(water.size());
      for (int k = 0; k < oxygen.nop; ++k) {
        oxygen.pts[static_cast<std::size_t>(k)] = work.pts[static_cast<std::size_t>(water[static_cast<std::size_t>(k)])];
        oxygen.pts[static_cast<std::size_t>(k)].type = 1;
        oxygen.pts[static_cast<std::size_t>(k)].atomID = k + 1;
        oxygen.idIndexMap[k + 1] = k;
      }
      if (!chill::chillPlusPublishedCutoff(oxygen, cutoff, 1)) {
        std::vector<int> offsets;
        std::vector<int> cols;
        nneigh::orthoFullCsr(oxygen, nullptr, cutoff, offsets, cols, false);
        chill::chillPlusFromCsr(oxygen, offsets, cols);
      }
      for (int k = 0; k < oxygen.nop; ++k) {
        labels[static_cast<std::size_t>(water[static_cast<std::size_t>(k)])] =
            static_cast<int>(oxygen.pts[static_cast<std::size_t>(k)].iceType);
      }
    }
    std::vector<bool> ice(static_cast<std::size_t>(work.nop), false);
    for (int index : water) {
      ice[static_cast<std::size_t>(index)] = crystal(labels[static_cast<std::size_t>(index)]);
    }
    site::FrontOptions options;
    const auto snap = site::analyseFront(work, ice, water, ions, options);
    ki("brine_found", snap.fronts.found ? 1 : 0);
    kd("brine_enter", snap.fronts.zEnter);
    kd("brine_exit", snap.fronts.zExit);
    kd("brine_thick", snap.fronts.thickness);
  }
  if (kind == "nucleation") {
    Cloud work = cloud;
    const auto rec = nucleation::analyzeFrame(work, cutoff, 7, 1, 0.36, 1);
    kd("nuc_mean", rec.q6barMean);
    kd("nuc_p50", rec.q6barP50);
    ki("nuc_ice", rec.q6barIce);
  }
  return 0;
}
