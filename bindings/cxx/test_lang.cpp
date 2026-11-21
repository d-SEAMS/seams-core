#include "seams_lang.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#if defined(_WIN32)
#include <windows.h>
#else
#include <dlfcn.h>
#endif

namespace {

using ArmFn = void (*)(int, int);
using LiveFn = long (*)();
using ResetFn = void (*)();

void *process_symbol(const char *name) {
#if defined(_WIN32)
  return reinterpret_cast<void *>(GetProcAddress(GetModuleHandleA(nullptr), name));
#else
  return dlsym(RTLD_DEFAULT, name);
#endif
}
ArmFn fault_arm() {
  return reinterpret_cast<ArmFn>(process_symbol("seams_test_arm"));
}
LiveFn fault_live() {
  return reinterpret_cast<LiveFn>(process_symbol("seams_test_live"));
}
LiveFn fault_peak() {
  return reinterpret_cast<LiveFn>(process_symbol("seams_test_peak"));
}
ResetFn fault_reset() {
  return reinterpret_cast<ResetFn>(process_symbol("seams_test_reset_peak"));
}

int fail(const std::string &message) {
  std::cerr << message;
  const char *err = seams_last_error();
  if (err != nullptr && err[0] != '\0') {
    std::cerr << ": " << err;
  }
  std::cerr << "\n";
  return 1;
}

int count_of(const int counts[9], int label) { return counts[label]; }

} // namespace

int main(int argc, char **argv) {
  const std::string root = argc > 1 ? argv[1] : ".";
  const std::string tip4p = root + "/input/traj/exampleTraj.lammpstrj";
  const std::string cubic = root + "/input/traj/mW_cubic.lammpstrj";
  const std::string si = root + "/input/traj/genice_sI.lammpstrj";

  if (std::string(seams_lang_version()).empty()) {
    return fail("empty version");
  }

  seams_frame *tip = seams_frame_read(tip4p.c_str(), 1, 2);
  if (tip == nullptr) {
    return fail("read tip4p");
  }
  if (seams_frame_nop(tip) != 250) {
    return fail("tip4p nop");
  }
  const double *xyz = seams_frame_positions(tip);
  if (xyz == nullptr) {
    return fail("positions");
  }

  std::vector<int> labels(250);
  int counts[9] = {};
  if (seams_frame_chill_plus(tip, 3.5, 2, "mutual", labels.data(), counts) !=
      0) {
    return fail("chill+ mutual");
  }
  if (count_of(counts, 4) != 1 || count_of(counts, 5) != 7 ||
      count_of(counts, 2) != 242) {
    std::cerr << "mutual chill+ " << counts[4] << " " << counts[5] << " "
              << counts[2] << "\n";
    return fail("mutual chill+ counts");
  }

  if (seams_frame_chill_plus(tip, 3.5, 2, "cutoff", labels.data(), counts) !=
      0) {
    return fail("chill+ cutoff");
  }
  if (count_of(counts, 5) != 9 || count_of(counts, 2) != 241 ||
      count_of(counts, 4) != 0) {
    std::cerr << "cutoff chill+ clathrate " << counts[4] << " inter "
              << counts[5] << " water " << counts[2] << "\n";
    return fail("cutoff chill+ counts");
  }
  int default_counts[9] = {};
  if (seams_frame_chill_plus(tip, 3.5, 2, nullptr, labels.data(),
                             default_counts) != 0) {
    return fail("chill+ default");
  }
  for (int i = 0; i < 9; ++i) {
    if (default_counts[i] != counts[i]) {
      return fail("chill+ default is not the cutoff list");
    }
  }
  if (seams_frame_chill_plus(tip, 3.5, 2, "knn-union", labels.data(),
                             counts) == 0) {
    return fail("chill+ accepted knn-union");
  }

  int n_hex = -1, n_cubic = -1, n_water = -1, n_both = -1;
  if (seams_frame_cages(tip, 3.5, 2, 4, "seeded", 0, nullptr, &n_hex, &n_cubic,
                        &n_water, &n_both) != 0) {
    return fail("cages");
  }
  if (n_hex != 0 || n_cubic != 0 || n_water != 250) {
    std::cerr << "cages " << n_hex << " " << n_cubic << " " << n_water << "\n";
    return fail("cage counts");
  }

  std::vector<double> ql(250), qlbar(250);
  if (seams_frame_steinhardt(tip, 3.5, 2, 3, ql.data(), qlbar.data()) != 0) {
    return fail("steinhardt");
  }
  double mean = 0.0;
  int nfinite = 0;
  for (double value : ql) {
    if (std::isfinite(value)) {
      mean += value;
      ++nfinite;
    }
  }
  mean = nfinite > 0 ? mean / static_cast<double>(nfinite) : 0.0;
  // Orthorhombic frames take the cutoff shell. 0.389362 is the
  // four-nearest mean of the same oxygens.
  if (std::abs(mean - 0.374897) > 1e-4) {
    std::cerr << "q3 mean " << mean << "\n";
    return fail("q3");
  }

  const int bins = 4;
  std::vector<double> r(bins), g(bins);
  std::vector<int> hist(bins);
  double volume = 0.0;
  if (seams_frame_rdf(tip, 6.0, bins, 2, 2, r.data(), g.data(), hist.data(),
                      &volume) != 0) {
    return fail("rdf");
  }
  if (hist[1] != 390 || hist[2] != 705 || hist[3] != 531) {
    std::cerr << "rdf counts " << hist[0] << " " << hist[1] << " " << hist[2]
              << " " << hist[3] << "\n";
    return fail("rdf counts");
  }

  std::vector<int> rings(9, 0);
  if (seams_frame_rings(tip, 3.5, 2, "mutual", 8, rings.data(), 9) != 0) {
    return fail("rings");
  }
  std::cout << "tip4p mutual rings";
  for (std::size_t i = 3; i <= 8; ++i) {
    std::cout << " " << i << ":" << rings[i];
  }
  std::cout << "\n";

  seams_frame *ice = seams_frame_read(cubic.c_str(), 1, 1);
  if (ice == nullptr || seams_frame_nop(ice) != 4096) {
    return fail("mW");
  }
  std::vector<int> ice_labels(4096);
  int ice_counts[9] = {};
  if (seams_frame_chill_plus(ice, 3.5, 1, "mutual", ice_labels.data(),
                             ice_counts) != 0) {
    return fail("mW chill+");
  }
  if (ice_counts[0] != 4096) {
    std::cerr << "mW cubic " << ice_counts[0] << "\n";
    return fail("mW cubic count");
  }
  if (seams_frame_cages(ice, 3.5, 1, 4, "seeded", 0, nullptr, &n_hex, &n_cubic,
                        &n_water, &n_both) != 0) {
    return fail("mW cages");
  }
  if (n_hex != 0 || n_cubic != 4096 || n_water != 0) {
    std::cerr << "mW cages " << n_hex << " " << n_cubic << " " << n_water
              << "\n";
    return fail("mW cage counts");
  }
  if (seams_path_nframes(cubic.c_str()) < 2) {
    return fail("mW nframes");
  }

  seams_frame *guest = seams_frame_read(si.c_str(), 1, 1);
  if (guest == nullptr) {
    return fail("sI");
  }
  int n_cages = 0;
  int n_atoms = 0;
  if (seams_frame_signature(guest, 3.5, 1, "cutoff", "512", &n_cages,
                            &n_atoms) != 0) {
    return fail("sI signature");
  }
  if (seams_frame_nop(guest) != 46 || n_cages != 2 || n_atoms != 40) {
    std::cerr << "sI " << seams_frame_nop(guest) << " cages " << n_cages
              << " atoms " << n_atoms << "\n";
    return fail("sI 512 census");
  }
  if (seams_frame_signature_edge(guest, 3.5, 1, "cutoff", "512", "none",
                                 &n_cages, &n_atoms) != 0 ||
      n_cages != 2 || n_atoms != 40) {
    std::cerr << "sI none cages " << n_cages << " atoms " << n_atoms << "\n";
    return fail("sI 512 none");
  }
  if (seams_frame_signature_edge(guest, 3.5, 1, "cutoff", "512", "1,2",
                                 &n_cages, &n_atoms) != 0 ||
      n_cages != 0) {
    std::cerr << "sI tight cages " << n_cages << "\n";
    return fail("sI 512 tight window");
  }
  if (seams_frame_signature(guest, 3.5, 1, "cutoff", "51268", &n_cages,
                            &n_atoms) != 0) {
    return fail("sI 51268");
  }
  if (n_cages != 0) {
    std::cerr << "sI 51268 cages " << n_cages << "\n";
    return fail("sI 51268 census");
  }

  const std::string nucleation = root + "/input/traj/nucleation_f1.lammpstrj";
  seams_frame *molten = seams_frame_read(nucleation.c_str(), 1, 1);
  if (molten == nullptr) {
    return fail("nucleation frame");
  }
  const int nMolten = seams_frame_nop(molten);
  std::vector<int> hostLabels(static_cast<std::size_t>(nMolten));
  std::vector<int> offLabels(static_cast<std::size_t>(nMolten));
  int hHex = 0, hCub = 0, hWat = 0, hBoth = 0;
  int oHex = 0, oCub = 0, oWat = 0, oBoth = 0;
  unsetenv("SEAMS_OFFLOAD");
  if (seams_frame_cages(molten, 3.5, 1, 4, "seeded", 0, hostLabels.data(),
                        &hHex, &hCub, &hWat, &hBoth) != 0) {
    return fail("seeded host");
  }
  setenv("SEAMS_OFFLOAD", "1", 1);
  if (seams_frame_cages(molten, 3.5, 1, 4, "seeded", 0, offLabels.data(),
                        &oHex, &oCub, &oWat, &oBoth) != 0) {
    return fail("seeded offload");
  }
  unsetenv("SEAMS_OFFLOAD");
  if (hHex != oHex || hCub != oCub || hWat != oWat || hBoth != oBoth ||
      hostLabels != offLabels) {
    std::cerr << "offload " << oHex << " " << oCub << " " << oBoth << " "
              << oWat << " host " << hHex << " " << hCub << " " << hBoth << " "
              << hWat << "\n";
    return fail("seeded offload disagrees with the host rule");
  }
  if (hHex + hCub + hBoth + hWat != nMolten) {
    std::cerr << "counters " << hHex << " " << hCub << " " << hBoth << " "
              << hWat << " nop " << nMolten << "\n";
    return fail("disjoint cage counters");
  }
  int codes[4] = {};
  for (int code : hostLabels) {
    if (code < 0 || code > 3) {
      return fail("cage code");
    }
    ++codes[code];
  }
  // 0 water, 1 hexagonal only, 2 cubic only, 3 both
  if (codes[1] != hHex || codes[2] != hCub || codes[3] != hBoth ||
      codes[0] != hWat || hBoth <= 0) {
    std::cerr << "codes " << codes[1] << " " << codes[2] << " " << codes[3]
              << " both " << hBoth << "\n";
    return fail("stacking-fault counters");
  }
  if (argc > 2) {
    const std::string cmd = std::string(argv[2]) + " cages " + nucleation +
                            " --graph seeded --cutoff 3.5";
#if defined(_WIN32)
    FILE *pipe = _popen(cmd.c_str(), "r");
#else
    FILE *pipe = popen(cmd.c_str(), "r");
#endif
    if (pipe == nullptr) {
      return fail("cli pipe");
    }
    std::string line;
    char buf[512];
    while (fgets(buf, sizeof buf, pipe) != nullptr) {
      line += buf;
    }
#if defined(_WIN32)
    const int rc = _pclose(pipe);
#else
    const int rc = pclose(pipe);
#endif
    if (rc != 0) {
      std::cerr << line;
      return fail("cli cages");
    }
    int cliNop = 0, cliHex = 0, cliCub = 0, cliBoth = 0, cliWat = 0;
    std::istringstream in(line);
    std::string tok;
    in >> tok >> cliNop >> tok >> tok >> tok >> cliHex >> tok >> cliCub >>
        tok >> cliBoth >> tok >> cliWat;
    if (cliHex != hHex || cliCub != hCub || cliBoth != hBoth ||
        cliWat != hWat) {
      std::cerr << "cli " << line << "abi " << hHex << " " << hCub << " "
                << hBoth << " " << hWat << "\n";
      return fail("abi disagrees with the cli");
    }
  }

  double box9[9] = {};
  seams_frame_box(tip, box9);
  if (!(box9[3] > 0.0)) {
    return fail("box");
  }
  if (seams_frame_index(tip) < 1) {
    return fail("index");
  }
  if (seams_frame_types(tip) == nullptr || seams_frame_ids(tip) == nullptr) {
    return fail("arrays");
  }
  if (seams_frame_positions(nullptr) != nullptr) {
    return fail("null positions");
  }
  if (seams_frame_chill(tip, 3.5, 2, "cutoff", labels.data(), counts) != 0) {
    return fail("plain chill");
  }
  if (seams_frame_chill_plus_rule(tip, 3.5, 2, "cutoff", "ovito-chill+", labels.data(), counts) != 0) {
    return fail("ovito rule");
  }
  if (seams_frame_chill_plus_rule(tip, 3.5, 2, "cutoff", "not-a-rule", labels.data(), counts) == 0) {
    return fail("bad rule accepted");
  }
  if (seams_frame_from_arrays(nullptr, 1, box9, nullptr, nullptr) != nullptr) {
    return fail("null xyz");
  }
  seams_frame *built = seams_frame_from_arrays(xyz, 4, box9, nullptr, nullptr);
  if (built == nullptr || seams_frame_nop(built) != 4) {
    return fail("from arrays");
  }
  seams_frame_free(built);

  const int graphOff[7] = {0, 1, 2, 3, 4, 5, 6};
  const int graphCol[6] = {1, 2, 3, 4, 5, 0};
  int graphHist[7] = {};
  int nGraphRings = 0;
  int nGraphVertices = 0;
  if (seams_graph_rings(graphOff, 7, graphCol, 6, 6, graphHist, 7,
                        nullptr, 0, nullptr, 0, &nGraphRings,
                        &nGraphVertices) != 0 ||
      graphHist[6] != 1 || nGraphRings != 1 || nGraphVertices != 6) {
    return fail("explicit graph ring query");
  }
  int packedOff[2] = {};
  int packedVertices[6] = {};
  if (seams_graph_rings(graphOff, 7, graphCol, 6, 6, graphHist, 7,
                        packedOff, 2, packedVertices, 6, &nGraphRings,
                        &nGraphVertices) != 0 ||
      packedOff[1] != 6) {
    return fail("explicit graph rings");
  }

  const double tetraXyz[18] = {
      50, 50, 50, 51, 51, 51, 51, 49, 49,
      49, 51, 49, 49, 49, 51, 50, 50, 52,
  };
  const double analyticBox[9] = {0, 0, 0, 100, 100, 100, 0, 0, 0};
  seams_frame *analytic =
      seams_frame_from_arrays(tetraXyz, 6, analyticBox, nullptr, nullptr);
  const int centre[1] = {0};
  const int tetraPartners[4] = {1, 2, 3, 4};
  double tetraQ[1] = {};
  if (analytic == nullptr ||
      seams_frame_tetrahedral(analytic, centre, 1, tetraPartners, 4, 1,
                              2.0, tetraQ) != 0 ||
      std::abs(tetraQ[0] - 1.0) > 1e-12) {
    return fail("language tetrahedral order");
  }
  seams_frame_free(analytic);

  const double radialXyz[15] = {
      50, 50, 50, 51, 50, 50, 50, 51.5, 50,
      50, 50, 52.5, 47, 50, 50,
  };
  const int radialPartners[4] = {1, 2, 3, 4};
  analytic =
      seams_frame_from_arrays(radialXyz, 5, analyticBox, nullptr, nullptr);
  double lsi[1] = {};
  if (analytic == nullptr ||
      seams_frame_lsi(analytic, centre, 1, radialPartners, 4, 3.1, lsi) !=
          0 ||
      std::abs(lsi[0] - 1.0 / 18.0) > 1e-12) {
    return fail("language LSI");
  }
  seams_frame_free(analytic);

  const double windingXyz[9] = {1, 4, 4, 4, 4, 4, 7, 4, 4};
  const double windingBox[9] = {0, 0, 0, 9, 9, 9, 0, 0, 0};
  const int triangleOff[4] = {0, 1, 2, 3};
  const int triangleCol[3] = {1, 2, 0};
  int triangleHist[4] = {};
  int trianglePackedOff[2] = {};
  int triangleVertices[3] = {};
  int triangleRings = 0;
  int triangleVertexCount = 0;
  analytic =
      seams_frame_from_arrays(windingXyz, 3, windingBox, nullptr, nullptr);
  if (analytic == nullptr ||
      seams_graph_rings(triangleOff, 4, triangleCol, 3, 3, triangleHist,
                        4, trianglePackedOff, 2, triangleVertices, 3,
                        &triangleRings, &triangleVertexCount) != 0) {
    return fail("language winding ring");
  }
  int axisWinding[1] = {};
  int isContractible[1] = {1};
  if (seams_frame_ring_winding(analytic, trianglePackedOff, triangleRings,
                               triangleVertices, triangleVertexCount, 0,
                               axisWinding, isContractible) != 0 ||
      axisWinding[0] == 0 || isContractible[0] != 0) {
    return fail("language ring winding");
  }
  seams_frame_free(analytic);

  if (seams_frame_prisms(nullptr, 3.5, 2, 6, nullptr, nullptr, nullptr) == 0) {
    return fail("null prisms");
  }
  int n_rows = 0;
  if (seams_frame_prisms(tip, 3.5, 2, 6, nullptr, nullptr, &n_rows) != 0 ||
      n_rows != 4) return fail("prism count");
  std::vector<seams_prism_row> prism_rows(static_cast<std::size_t>(n_rows));
  if (seams_frame_prisms(tip, 3.5, 2, 6, nullptr, prism_rows.data(), &n_rows) != 0) {
    return fail("prisms");
  }
  const int nbin = seams_frame_rdf2d_nbin(3.5, 0.0);
  std::vector<double> r2(static_cast<std::size_t>(nbin));
  std::vector<double> g2(static_cast<std::size_t>(nbin));
  if (seams_frame_rdf2d(nullptr, 3.5, 0.5, nullptr, r2.data(), g2.data()) == 0) {
    return fail("null rdf2d");
  }
  if (seams_frame_rdf2d(tip, 3.5, 0.5, nullptr, r2.data(), g2.data()) != 0) {
    return fail("rdf2d");
  }
  int sheet[8] = {};
  if (seams_frame_monolayer(tip, 3.5, 2, 6, 1.0, nullptr, sheet, 8) != 0) {
    return fail("monolayer");
  }
  if (seams_frame_chill(nullptr, 3.5, 2, "cutoff", labels.data(), counts) == 0) {
    return fail("null chill");
  }
  if (seams_last_error() == nullptr) {
    return fail("error string");
  }

  if (seams_path_nframes(nullptr) != -1) {
    return fail("null path");
  }
  if (seams_path_nframes("missing.xyz") != -1) {
    return fail("xyz frames");
  }
  if (std::string(seams_last_error()).find("cannot read") == std::string::npos) {
    return fail("xyz text");
  }
  if (seams_frame_read(nullptr, 1, 1) != nullptr) {
    return fail("null read");
  }
  seams_frame *missing = seams_frame_read("missing-lang.xyz", 1, 1);
  if (missing != nullptr) {
    return fail("missing xyz opened");
  }
  missing = seams_frame_read("missing-lang.pdb", 1, 1);
  if (missing != nullptr) {
    seams_frame_free(missing);
  }
  missing = seams_frame_read("missing-lang.gro", 1, 0);
  if (missing != nullptr) {
    seams_frame_free(missing);
  }
  missing = seams_frame_read("missing-lang.dcd", 1, -1);
  if (missing != nullptr) {
    seams_frame_free(missing);
  }
  missing = seams_frame_read("missing-lang.con", 1, 1);
  if (missing != nullptr) {
    seams_frame_free(missing);
  }
  seams_frame *allTypes = seams_frame_read(tip4p.c_str(), 1, -1);
  if (allTypes == nullptr) {
    return fail("type filter");
  }
  int hb_n = -1;
  if (seams_frame_hbond_edges(allTypes, 3.5, 2, 1, 2.42, 30.0, 0, nullptr, nullptr, nullptr,
                              nullptr, nullptr, 0, &hb_n) != 0) {
    seams_frame_free(allTypes);
    return fail("hbond edge query");
  }
  if (hb_n < 0) {
    seams_frame_free(allTypes);
    return fail("hbond edge count");
  }
  std::vector<uint32_t> donor(static_cast<std::size_t>(std::max(hb_n, 1)));
  std::vector<uint32_t> acceptor(static_cast<std::size_t>(std::max(hb_n, 1)));
  std::vector<uint32_t> hydrogen(static_cast<std::size_t>(std::max(hb_n, 1)));
  std::vector<double> dDA(static_cast<std::size_t>(std::max(hb_n, 1)));
  std::vector<double> angleDeg(static_cast<std::size_t>(std::max(hb_n, 1)));
  int hb_n_copy = hb_n;
  if (seams_frame_hbond_edges(allTypes, 3.5, 2, 1, 2.42, 30.0, 3, donor.data(), acceptor.data(),
                              hydrogen.data(), dDA.data(), angleDeg.data(), hb_n_copy,
                              &hb_n_copy) != 0) {
    seams_frame_free(allTypes);
    return fail("hbond edge copy");
  }
  if (hb_n_copy != hb_n) {
    seams_frame_free(allTypes);
    return fail("hbond edge count changed");
  }
  if (hb_n > 0) {
    int hb_short = hb_n;
    if (seams_frame_hbond_edges(allTypes, 3.5, 2, 1, 2.42, 30.0, 3, donor.data(),
                                acceptor.data(), hydrogen.data(), dDA.data(), angleDeg.data(),
                                hb_n - 1, &hb_short) == 0) {
      seams_frame_free(allTypes);
      return fail("hbond edge short buffer accepted");
    }
  }
  for (int i = 1; i < hb_n; ++i) {
    const bool sorted = (donor[static_cast<std::size_t>(i - 1)] <
                         donor[static_cast<std::size_t>(i)]) ||
                        (donor[static_cast<std::size_t>(i - 1)] ==
                             donor[static_cast<std::size_t>(i)] &&
                         acceptor[static_cast<std::size_t>(i - 1)] <
                             acceptor[static_cast<std::size_t>(i)]) ||
                        (donor[static_cast<std::size_t>(i - 1)] ==
                             donor[static_cast<std::size_t>(i)] &&
                         acceptor[static_cast<std::size_t>(i - 1)] ==
                             acceptor[static_cast<std::size_t>(i)] &&
                         hydrogen[static_cast<std::size_t>(i - 1)] <=
                             hydrogen[static_cast<std::size_t>(i)]);
    if (!sorted) {
      seams_frame_free(allTypes);
      return fail("hbond edge order");
    }
  }
  seams_frame_free(allTypes);
  seams_frame *guess = seams_frame_read(tip4p.c_str(), 1, 0);
  if (guess == nullptr) {
    return fail("type guess");
  }
  seams_frame_free(guess);

  {
    const double mix_xyz[18] = {
        0.0,  0.0,  0.0,  0.96, 0.0,  0.0, 0.0, 0.96, 0.0,
        2.8,  0.0,  0.0,  3.76, 0.0,  0.0, 2.8, 0.96, 0.0,
    };
    const double mix_box[9] = {0.0, 0.0, 0.0, 20.0, 20.0, 20.0, 0.0, 0.0, 0.0};
    const int mix_types[6] = {2, 1, 1, 3, 1, 1};
    const int mix_mols[6] = {1, 1, 1, 2, 2, 2};
    const int mix_ids[6] = {1, 2, 3, 4, 5, 6};
    seams_frame *mix_heavy =
        seams_frame_from_species(mix_xyz, 6, mix_box, mix_types, mix_ids, mix_mols);
    if (mix_heavy == nullptr) {
      return fail("mixed heavy frame");
    }
    int mix_n = -1;
    if (seams_frame_hbond_edges(mix_heavy, 3.5, 0, 1, 2.42, 30.0, 3, nullptr, nullptr, nullptr,
                                nullptr, nullptr, 0, &mix_n) != 0) {
      seams_frame_free(mix_heavy);
      return fail("mixed heavy hbond query");
    }
    if (mix_n != 1) {
      seams_frame_free(mix_heavy);
      return fail("mixed heavy hbond count");
    }
    uint32_t md = 0, ma = 0, mh_idx = 0;
    double mda = 0.0, mang = 0.0;
    if (seams_frame_hbond_edges(mix_heavy, 3.5, 0, 1, 2.42, 30.0, 3, &md, &ma, &mh_idx, &mda,
                                &mang, 1, &mix_n) != 0 ||
        md != 0 || ma != 3 || mh_idx != 1) {
      seams_frame_free(mix_heavy);
      return fail("mixed heavy hbond edge indices");
    }
    seams_frame_free(mix_heavy);
  }

  int types4[4] = {2, 2, 1, 1};
  int ids4[4] = {10, 11, 12, 13};
  double sheared[9] = {0, 0, 0, 10, 10, 10, 1.5, 0, 0};
  if (seams_frame_from_arrays(xyz, 0, box9, nullptr, nullptr) != nullptr) {
    return fail("zero count");
  }
  if (seams_frame_from_arrays(xyz, 4, nullptr, nullptr, nullptr) != nullptr) {
    return fail("null box");
  }
  seams_frame *typed = seams_frame_from_arrays(xyz, 4, sheared, types4, ids4);
  if (typed == nullptr || seams_frame_nop(typed) != 4) {
    return fail("typed arrays");
  }
  double filled[9] = {};
  seams_frame_box(typed, filled);
  if (!(filled[6] != 0.0)) {
    return fail("shear");
  }
  seams_frame_free(typed);
  seams_frame_free(nullptr);
  if (seams_frame_nop(nullptr) != 0) {
    return fail("null nop");
  }
  if (seams_frame_index(nullptr) != 0) {
    return fail("null index");
  }
  seams_frame_box(nullptr, box9);
  if (seams_frame_types(nullptr) != nullptr) {
    return fail("null types");
  }
  if (seams_frame_ids(nullptr) != nullptr) {
    return fail("null ids");
  }

  if (seams_frame_chill_plus(tip, 3.5, 2, "", labels.data(), counts) != 0) {
    return fail("empty graph");
  }
  if (seams_frame_chill_plus(tip, 3.5, 2, "cutoff", labels.data(), nullptr) != 0) {
    return fail("null counts");
  }
  if (seams_frame_chill_plus(tip, 3.5, 2, "knn", labels.data(), counts) == 0) {
    return fail("knn accepted");
  }
  std::string longRule(300, 'q');
  if (seams_frame_chill_plus_rule(tip, 3.5, 2, "cutoff", longRule.c_str(), labels.data(), counts) == 0) {
    return fail("long rule accepted");
  }
  if (seams_frame_rings(tip, 3.5, 2, "cutoff", 6, rings.data(), 9) != 0) {
    return fail("cutoff rings");
  }
  if (seams_frame_rings(tip, 3.5, 2, "knn-union", 6, rings.data(), 9) != 0) {
    return fail("union rings");
  }
  if (seams_frame_rings(nullptr, 3.5, 2, "mutual", 6, rings.data(), 9) == 0) {
    return fail("null rings");
  }
  if (seams_frame_rings(tip, 3.5, 2, "mutual", 2, rings.data(), 9) == 0) {
    return fail("shallow rings");
  }

  int cageHex = 0;
  if (seams_frame_cages(tip, 3.5, 2, 4, "cutoff", 0, nullptr, &cageHex, nullptr, nullptr, nullptr) != 0) {
    return fail("cutoff cages");
  }
  if (seams_frame_cages(nullptr, 3.5, 2, 4, "cutoff", 0, nullptr, &cageHex, nullptr, nullptr, nullptr) == 0) {
    return fail("null cages");
  }
  if (seams_frame_cages(tip, 3.5, 2, 4, "banana", 0, nullptr, &cageHex, nullptr, nullptr, nullptr) == 0) {
    return fail("banana cages");
  }
  if (seams_frame_cages(tip, 3.5, 2, 0, "seeded", 0, nullptr, &cageHex, nullptr, nullptr, nullptr) != 0) {
    return fail("k zero");
  }
  setenv("SEAMS_OFFLOAD", "1", 1);
  if (seams_frame_cages(tip, 3.5, 2, 4, "mutual", 0, nullptr, &cageHex, &cageHex, &cageHex, &cageHex) != 0) {
    return fail("offload cages");
  }
  unsetenv("SEAMS_OFFLOAD");
  if (seams_frame_steinhardt(nullptr, 3.5, 2, 4, ql.data(), qlbar.data()) == 0) {
    return fail("null steinhardt");
  }
  if (seams_frame_steinhardt(tip, 3.5, 2, 4, nullptr, nullptr) != 0) {
    return fail("steinhardt counts");
  }
  if (seams_frame_signature(nullptr, 3.5, 1, "cutoff", "512", &n_cages, &n_atoms) == 0) {
    return fail("null signature");
  }
  if (seams_frame_signature(guest, 3.5, 1, nullptr, "512", &n_cages, &n_atoms) != 0) {
    return fail("default signature graph");
  }
  if (seams_frame_rdf_nbins(6.0, 0) < 1) {
    return fail("default bins");
  }
  if (seams_frame_rdf(nullptr, 6.0, 4, 2, 2, r.data(), g.data(), hist.data(), &volume) == 0) {
    return fail("null rdf");
  }
  if (seams_frame_rdf(tip, 6.0, 4, 2, 0, r.data(), g.data(), hist.data(), nullptr) != 0) {
    return fail("rdf no volume");
  }

  const auto out = std::filesystem::temp_directory_path() / "dseams-lang-out";
  std::filesystem::create_directories(out);
  const std::string outDir = out.string();
  if (seams_frame_rdf2d(tip, 3.5, 0.5, outDir.c_str(), r2.data(), g2.data()) != 0) {
    return fail("rdf2d file");
  }
  if (seams_frame_rdf2d_nbin(3.5, -1.0) < 1) {
    return fail("rdf2d width");
  }
  if (seams_frame_prisms(tip, 3.5, 2, 6, outDir.c_str(), prism_rows.data(), &n_rows) != 0) {
    return fail("prism file");
  }
  if (seams_frame_prisms(tip, 3.5, 2, 6, nullptr, prism_rows.data(), nullptr) == 0) {
    return fail("short prism buffer");
  }
  int sheetWide[12] = {};
  if (seams_frame_monolayer(tip, 3.5, 2, 6, 1.0, outDir.c_str(), sheetWide, 12) != 0) {
    return fail("monolayer file");
  }
  if (seams_frame_monolayer(nullptr, 3.5, 2, 6, 1.0, nullptr, sheetWide, 12) == 0) {
    return fail("null monolayer");
  }
  if (seams_frame_monolayer(tip, 3.5, 2, 6, 1.0, nullptr, sheetWide, 2) == 0) {
    return fail("short monolayer");
  }
  const std::string ice1h = root + "/input/traj/ice1h.lammpstrj";
  seams_frame *none = seams_frame_read(ice1h.c_str(), 1, 9);
  if (none != nullptr) {
    if (seams_frame_nop(none) == 0 && seams_frame_chill_plus(none, 3.5, 0, "cutoff", nullptr, nullptr) != 0) {
      return fail("filtered empty");
    }
    seams_frame_free(none);
  }
  if (seams_frame_chill_plus(tip, 3.5, 1, "cutoff", labels.data(), counts) != 0) {
    return fail("other type");
  }
  seams_frame *guessed = seams_frame_read(ice1h.c_str(), 1, 0);
  if (guessed == nullptr || seams_frame_nop(guessed) == 0) {
    return fail("ice guess");
  }
  seams_frame_free(guessed);
  const auto emptyPath = out / "empty.lammpstrj";
  {
    std::ofstream empty(emptyPath);
    empty << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n0\n"
             "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
             "ITEM: ATOMS id type x y z\n";
  }
  seams_frame *blank = seams_frame_read(emptyPath.c_str(), 1, 0);
  if (blank == nullptr) {
    return fail("empty frame");
  }
  if (seams_path_nframes(emptyPath.c_str()) != 1 ||
      seams_frame_read(emptyPath.c_str(), 2, 0) != nullptr ||
      seams_path_nframes(emptyPath.c_str()) != 1) {
    seams_frame_free(blank);
    return fail("past-end frame read");
  }
  if (seams_frame_chill_plus(blank, 3.5, 0, "cutoff", nullptr, nullptr) != 0) {
    return fail("empty chill");
  }
  if (seams_frame_cages(blank, 3.5, 0, 4, "cutoff", 0, nullptr, &cageHex, nullptr, nullptr, nullptr) != 0) {
    return fail("empty cages");
  }
  seams_frame_free(blank);
  if (seams_frame_rings(tip, 3.5, 2, "", 6, rings.data(), 9) != 0) {
    return fail("empty-name rings");
  }
  if (seams_frame_rings(tip, 3.5, 2, "banana", 6, rings.data(), 9) == 0) {
    return fail("banana rings");
  }
  if (seams_frame_rings(tip, 3.5, 2, "seeded", 6, rings.data(), 9) == 0) {
    return fail("seeded rings");
  }
  if (seams_frame_signature(guest, 3.5, 1, "cutoff", "not-a-signature", &n_cages, &n_atoms) == 0) {
    return fail("bad signature");
  }
  const auto corrupt = out / "corrupt.lammpstrj";
  {
    std::ofstream bad(corrupt);
    bad << "not a lammpstrj\n";
  }
  if (seams_path_nframes(corrupt.c_str()) > 0) {
    return fail("corrupt frames");
  }
  if (seams_path_nframes("no-extension") > 0) {
    return fail("bare frames");
  }
  const std::string xyzPath = root + "/input/xyz/configuration.xyz";
  if (seams_path_nframes(xyzPath.c_str()) != 1) {
    return fail("xyz frames");
  }
  seams_frame *xyzFrame = seams_frame_read(xyzPath.c_str(), 1, 1);
  if (xyzFrame != nullptr) {
    double xyzBox[9] = {};
    seams_frame_box(xyzFrame, xyzBox);
    seams_frame_free(xyzFrame);
  }
  if (seams_frame_rings(tip, 3.5, 2, "knn-mutual", 6, rings.data(), 9) != 0) {
    return fail("knn rings");
  }
  if (seams_frame_rings(tip, 3.5, 2, "union", 6, rings.data(), 9) != 0) {
    return fail("union name");
  }
  if (seams_frame_rings(tip, 3.5, 2, nullptr, 6, rings.data(), 9) != 0) {
    return fail("null graph rings");
  }
  double shearY[9] = {0, 0, 0, 10, 10, 10, 0, 1.0, 0};
  double shearZ[9] = {0, 0, 0, 10, 10, 10, 0, 0, 1.0};
  seams_frame *sy = seams_frame_from_arrays(xyz, 4, shearY, types4, ids4);
  seams_frame *sz = seams_frame_from_arrays(xyz, 4, shearZ, nullptr, nullptr);
  if (sy == nullptr || sz == nullptr) {
    return fail("shear axes");
  }
  seams_frame_free(sy);
  seams_frame_free(sz);
  std::vector<int> per(static_cast<std::size_t>(seams_frame_nop(tip)));
  if (seams_frame_cages(tip, 3.5, 2, 4, nullptr, 1, per.data(), nullptr, nullptr, nullptr, nullptr) != 0) {
    return fail("null graph cages");
  }
  if (seams_frame_signature(guest, 3.5, 1, "", "512", &n_cages, nullptr) == 0) {
    return fail("null atom count");
  }
  if (seams_frame_rdf(tip, 6.0, 4, 2, 2, nullptr, g.data(), hist.data(), &volume) == 0) {
    return fail("null rdf r");
  }
  const std::string dir = (out / "sub").string();
  std::filesystem::create_directories(dir);
  if (seams_frame_rdf2d(tip, 3.5, 0.5, dir.c_str(), r2.data(), g2.data()) != 0) {
    return fail("rdf2d slash");
  }
  if (seams_frame_prisms(tip, 3.5, 2, 4, dir.c_str(), nullptr, &n_rows) != 0) {
    return fail("prism dir");
  }
  if (seams_frame_monolayer(tip, 3.5, 2, 4, 1.0, dir.c_str(), sheetWide, 12) != 0) {
    return fail("sheet dir");
  }

  seams_frame *oxygen = seams_frame_read(ice1h.c_str(), 1, 0);
  if (oxygen == nullptr || seams_frame_nop(oxygen) == 0) {
    return fail("oxygen guess");
  }
  if (seams_frame_chill(oxygen, 3.5, 0, "cutoff", labels.data(), counts) != 0) {
    return fail("chill type 0");
  }
  if (seams_frame_chill(oxygen, 3.5, 1, nullptr, nullptr, nullptr) != 0) {
    return fail("chill null graph");
  }
  if (seams_frame_chill(nullptr, 3.5, 1, "cutoff", nullptr, nullptr) == 0) {
    return fail("chill null frame");
  }
  if (seams_frame_chill_plus(oxygen, 3.5, 2, "mutual", labels.data(), counts) != 0) {
    return fail("skip chill");
  }
  if (seams_frame_positions(nullptr) != nullptr) {
    return fail("null positions");
  }
  if (seams_frame_types(nullptr) != nullptr) {
    return fail("null types");
  }
  if (seams_frame_ids(nullptr) != nullptr) {
    return fail("null ids");
  }
  if (seams_frame_positions(oxygen) == nullptr) {
    return fail("oxygen positions");
  }
  if (seams_frame_types(oxygen) == nullptr) {
    return fail("oxygen types");
  }
  if (seams_frame_ids(oxygen) == nullptr) {
    return fail("oxygen ids");
  }
  double ignoredBox[9] = {};
  seams_frame_box(nullptr, ignoredBox);
  seams_frame_box(oxygen, nullptr);
  seams_frame *blank2 = seams_frame_read(emptyPath.c_str(), 1, -1);
  if (blank2 == nullptr) {
    return fail("blank again");
  }
  if (seams_frame_positions(blank2) != nullptr) {
    return fail("empty positions");
  }
  if (seams_frame_types(blank2) != nullptr) {
    return fail("empty types");
  }
  if (seams_frame_ids(blank2) != nullptr) {
    return fail("empty ids");
  }
  int emptyLabels[1] = {};
  if (seams_frame_chill_plus(blank2, 3.5, 0, "cutoff", emptyLabels, counts) != 0) {
    return fail("empty labels");
  }
  if (seams_frame_cages(blank2, 3.5, 0, 4, "bogus", 0, nullptr, &cageHex,
                        nullptr, nullptr, nullptr) == 0) {
    return fail("empty unknown cage graph");
  }
  if (seams_frame_rdf2d(blank2, 3.5, 0.5, nullptr, r2.data(), g2.data()) != 0) {
    return fail("empty rdf2d");
  }
  for (double value : g2) {
    if (value != 0.0) {
      return fail("empty rdf2d values");
    }
  }
  int emptyPrismN = 4;
  seams_prism_row emptyPrisms[4] = {};
  if (seams_frame_prisms(blank2, 3.5, 0, 6, outDir.c_str(), emptyPrisms,
                         &emptyPrismN) != 0) {
    return fail("empty prisms");
  }
  int perEmpty[1] = {};
  if (seams_frame_cages(blank2, 3.5, 0, 4, "cutoff", 0, perEmpty, &cageHex, nullptr, nullptr, nullptr) != 0) {
    return fail("empty per atom");
  }
  if (seams_frame_cages(oxygen, 3.5, 2, 4, "cutoff", 0, nullptr, &cageHex, nullptr, nullptr, nullptr) != 0) {
    return fail("skip type");
  }
  int wide[16] = {};
  if (seams_frame_rings(oxygen, 3.5, 1, "knn", 6, wide, 16) != 0) {
    return fail("knn rings");
  }
  int ncages = 0;
  int natoms = 0;
  if (seams_frame_signature(oxygen, 3.5, 1, "", "512", &ncages, &natoms) != 0) {
    return fail("empty graph signature");
  }
  if (seams_frame_signature_edge(oxygen, 3.5, 1, "cutoff", "512", "none", &ncages, &natoms) != 0) {
    return fail("edge none");
  }
  if (seams_frame_signature(oxygen, 3.5, 1, "mutual", "512", &ncages, &natoms) != 0) {
    return fail("mutual signature");
  }
  if (seams_frame_signature(nullptr, 3.5, 1, "cutoff", "512", &ncages, &natoms) == 0) {
    return fail("null signature frame");
  }
  if (seams_frame_signature(oxygen, 3.5, 1, "cutoff", nullptr, &ncages, &natoms) == 0) {
    return fail("null census");
  }
  double rr[8] = {};
  double gg[8] = {};
  int cc[8] = {};
  double vol = 0.0;
  if (seams_frame_rdf(oxygen, 6.0, 4, 0, 0, rr, gg, cc, &vol) != 0) {
    return fail("rdf type 0");
  }
  if (seams_frame_rdf(oxygen, 6.0, 4, 1, -1, rr, gg, cc, nullptr) != 0) {
    return fail("rdf type j");
  }
  if (seams_frame_rdf(oxygen, 6.0, 4, 1, 1, rr, nullptr, cc, &vol) == 0) {
    return fail("rdf null g");
  }
  if (seams_frame_rdf(oxygen, 6.0, 4, 1, 1, rr, gg, nullptr, &vol) == 0) {
    return fail("rdf null count");
  }
  if (seams_frame_rdf(blank2, 6.0, 0, 0, 0, rr, gg, cc, &vol) == 0) {
    return fail("rdf zero bins");
  }
  if (std::string(seams_last_error()).find("positive bin count") ==
      std::string::npos) {
    return fail("rdf zero bins");
  }
  if (seams_frame_rdf_nbins(6.0, 12) != 12) {
    return fail("rdf bin query");
  }
  const int autoBins = seams_frame_rdf_nbins(6.0, 0);
  std::vector<double> sizedR(static_cast<std::size_t>(autoBins));
  std::vector<double> sizedG(static_cast<std::size_t>(autoBins));
  std::vector<int> sizedC(static_cast<std::size_t>(autoBins));
  if (seams_frame_rdf(blank2, 6.0, autoBins, 0, 0, sizedR.data(), sizedG.data(),
                     sizedC.data(), &vol) != 0) {
    return fail("rdf sized bins");
  }
  constexpr int wideBins = 64;
  std::vector<double> wideR(wideBins);
  std::vector<double> wideG(wideBins);
  std::vector<int> wideC(wideBins);
  if (seams_frame_rdf(blank2, 6.0, wideBins, 0, 0, wideR.data(), wideG.data(),
                     wideC.data(), &vol) != 0) {
    return fail("rdf wide bins");
  }
  if (seams_frame_rdf2d(oxygen, 3.5, 0.5, nullptr, rr, gg) != 0) {
    return fail("rdf2d memory");
  }
  if (seams_frame_rdf2d(oxygen, 3.5, 0.5, "", rr, gg) != 0) {
    return fail("rdf2d empty dir");
  }
  const std::string slash = outDir + "/";
  if (seams_frame_rdf2d(oxygen, 3.5, 0.5, slash.c_str(), rr, gg) != 0) {
    return fail("rdf2d slash");
  }
  if (seams_frame_rdf2d(nullptr, 3.5, 0.5, nullptr, rr, gg) == 0) {
    return fail("rdf2d null frame");
  }
  if (seams_frame_rdf2d(oxygen, 3.5, 0.5, nullptr, nullptr, gg) == 0) {
    return fail("rdf2d null r");
  }
  if (seams_frame_rdf2d(oxygen, 3.5, 0.5, nullptr, rr, nullptr) == 0) {
    return fail("rdf2d null g");
  }
  int sheetCounts[16] = {};
  if (seams_frame_monolayer(oxygen, 3.5, 1, 6, 1.0, nullptr, sheetCounts, 16) != 0) {
    return fail("sheet memory");
  }
  if (seams_frame_monolayer(oxygen, 3.5, 1, 6, 1.0, "", sheetCounts, 16) != 0) {
    return fail("sheet empty dir");
  }
  if (seams_frame_monolayer(oxygen, 3.5, 1, 6, 1.0, slash.c_str(), sheetCounts, 16) != 0) {
    return fail("sheet slash");
  }
  std::vector<double> qlOxygen(static_cast<std::size_t>(seams_frame_nop(oxygen)));
  std::vector<double> qbarOxygen(qlOxygen.size());
  if (seams_frame_steinhardt(oxygen, 3.5, 0, 4, qlOxygen.data(), qbarOxygen.data()) != 0) {
    return fail("steinhardt type 0");
  }
  if (seams_frame_steinhardt(oxygen, 3.5, 1, 4, nullptr, qbarOxygen.data()) != 0) {
    return fail("steinhardt null ql");
  }
  if (seams_frame_steinhardt(oxygen, 3.5, 1, 4, qlOxygen.data(), nullptr) != 0) {
    return fail("steinhardt null bar");
  }
  if (seams_frame_steinhardt(oxygen, 3.5, 1, -1, qlOxygen.data(), qbarOxygen.data()) != 0 && seams_last_error()[0] == '\0') {
    return fail("silent bad order");
  }
  if (seams_frame_steinhardt(oxygen, 3.5, 1, 5, qlOxygen.data(),
                              qbarOxygen.data()) == 0) {
    return fail("unsupported Steinhardt degree");
  }
  if (seams_frame_rings(oxygen, std::nan(""), 1, "cutoff", 6, wide, 16) == 0 ||
      seams_frame_rings(oxygen, 0.0, 1, "cutoff", 6, wide, 16) == 0) {
    return fail("invalid ring cutoff");
  }
  if (seams_frame_steinhardt(blank2, 3.5, 0, 4, qlOxygen.data(), qbarOxygen.data()) != 0) {
    return fail("empty steinhardt");
  }
  const std::string missingDir = (out / "no-such" / "dir").string();
  if (seams_frame_rdf2d(oxygen, 3.5, 0.5, missingDir.c_str(), rr, gg) == 0 && seams_last_error()[0] != '\0') {
    return fail("rdf2d error with success");
  }
  int prismCount = 8;
  std::vector<seams_prism_row> prismWide(8);
  if (seams_frame_prisms(oxygen, 3.5, 1, 4, missingDir.c_str(), prismWide.data(), &prismCount) == 0 && seams_last_error()[0] != '\0') {
    return fail("prism error with success");
  }
  if (seams_frame_monolayer(oxygen, 3.5, 1, 4, 1.0, missingDir.c_str(), sheetCounts, 16) == 0 && seams_last_error()[0] != '\0') {
    return fail("sheet error with success");
  }
  prismCount = 8;
  if (seams_frame_prisms(oxygen, 3.5, 1, 6, nullptr, prismWide.data(), &prismCount) != 0) {
    return fail("prism memory");
  }
  int shortRows = 0;
  if (seams_frame_prisms(oxygen, 3.5, 1, 6, nullptr, prismWide.data(), &shortRows) == 0) {
    return fail("prism short rows");
  }
  if (seams_frame_prisms(oxygen, 3.5, 1, 4, nullptr, prismWide.data(), nullptr) == 0) {
    return fail("prism null count");
  }
  if (seams_frame_prisms(oxygen, 3.5, 1, 4, slash.c_str(), nullptr, &prismCount) != 0) {
    return fail("prism slash");
  }
  if (seams_path_nframes("missing-frames.lammpstrj") > 0) {
    return fail("missing frames");
  }
  if (seams_path_nframes("missing-frames.xyz") != -1) {
    return fail("xyz frames");
  }
  if (std::string(seams_last_error()).find("cannot read") == std::string::npos) {
    return fail("xyz error");
  }
  if (seams_frame_read("missing-frames.xyz", 1, 1) != nullptr) {
    return fail("missing xyz read");
  }
  if (seams_frame_read("missing-frames.con", 1, 1) != nullptr) {
    return fail("missing con read");
  }
  if (seams_last_error() == nullptr || seams_last_error()[0] == '\0') {
    return fail("missing con error");
  }
  double pairXyz[6] = {0, 0, 0, 1, 0, 0};
  double pairBox[9] = {0, 0, 0, 5, 5, 5, 0, 0, 0};
  int pairTypes[2] = {1, 2};
  seams_frame *pair = seams_frame_from_arrays(pairXyz, 2, pairBox, pairTypes, nullptr);
  if (pair == nullptr) {
    return fail("pair");
  }
  int pairRings[8] = {};
  if (seams_frame_rings(pair, 3.5, 0, "mutual", 6, pairRings, 8) != 0) {
    return fail("pair rings");
  }
  int pairPrismN = 8;
  std::vector<seams_prism_row> pairPrism(8);
  if (seams_frame_prisms(pair, 3.5, 0, 6, nullptr, pairPrism.data(), &pairPrismN) != 0) {
    return fail("pair prisms");
  }
  if (seams_frame_cages(pair, 3.5, 1, 4, "seeded", 1, nullptr, &cageHex, &cageHex, &cageHex, &cageHex) != 0) {
    return fail("pair seeded");
  }
  if (seams_frame_cages(pair, 3.5, 0, 4, "cutoff", 0, nullptr, &cageHex, nullptr, nullptr, nullptr) != 0) {
    return fail("pair cutoff cages");
  }
  if (seams_frame_monolayer(pair, 3.5, 0, 6, 1.0, nullptr, sheetCounts, 16) != 0) {
    return fail("pair sheet");
  }
  if (seams_frame_cages(oxygen, 3.5, 1, 4, "cutoff", 0, nullptr, &cageHex, &cageHex, &cageHex, &cageHex) != 0) {
    return fail("hex cages");
  }
  const std::string cubicPath = root + "/input/traj/mW_cubic.lammpstrj";
  seams_frame *cubicFrame = seams_frame_read(cubicPath.c_str(), 1, 1);
  if (cubicFrame == nullptr) {
    return fail("cubic read");
  }
  if (seams_frame_cages(cubicFrame, 3.5, 1, 4, "cutoff", 0, nullptr, &cageHex, &cageHex, &cageHex, &cageHex) != 0) {
    return fail("cubic cages");
  }
  seams_frame_free(cubicFrame);
  double faultXyz[3] = {0, 0, 0};
  double faultBox[9] = {0, 0, 0, 10, 10, 10, 0, 0, 0};
  int faultRings[8] = {};
  double faultR[4] = {};
  double faultG[4] = {};
  int faultC[4] = {};
  int faultSheet[8] = {};
  seams_prism_row faultPrism[4] = {};
  int faultPrismN = 4;
  int faultCages = 0;
  int faultAtoms = 0;
  // A sanitizer build does not preload this allocator. The checks
  // below need it. The ordinary build does preload it.
  if (fault_arm() != nullptr && fault_live() != nullptr) {
  auto caught = [](const char *label, auto &&call) -> int {
    const long before = fault_live()();
    fault_reset()();
    fault_arm()(1, 1);
    const int code = call();
    fault_arm()(0, 0);
    if (code == 0) {
      return fail(label);
    }
    if (std::string(seams_last_error()).find("bad_alloc") == std::string::npos) {
      return fail(label);
    }
    if (fault_live()() != before) {
      return fail(label);
    }
    if (fault_peak()() <= before) {
      return fail(label);
    }
    return 0;
  };
  const std::string longMissing(80, 'q');
  if (caught("fault frames", [&] {
        return seams_path_nframes((longMissing + ".lammpstrj").c_str()) == -1
                   ? 1
                   : 0;
      })) {
    return 1;
  }
  if (caught("fault read", [&] {
        return seams_frame_read(ice1h.c_str(), 1, 1) == nullptr ? 1 : 0;
      })) {
    return 1;
  }
  if (caught("fault arrays", [&] {
        return seams_frame_from_arrays(faultXyz, 1, faultBox, nullptr,
                                       nullptr) == nullptr
                   ? 1
                   : 0;
      })) {
    return 1;
  }
  if (caught("fault chill", [&] {
        return seams_frame_chill(oxygen, 3.5, 1, "cutoff", nullptr, nullptr) !=
                       0
                   ? 1
                   : 0;
      })) {
    return 1;
  }
  if (caught("fault rings", [&] {
        return seams_frame_rings(oxygen, 3.5, 1, "mutual", 6, faultRings, 8) !=
                       0
                   ? 1
                   : 0;
      })) {
    return 1;
  }
  if (caught("fault steinhardt", [&] {
        return seams_frame_steinhardt(oxygen, 3.5, 1, 4, faultR, faultG) != 0
                   ? 1
                   : 0;
      })) {
    return 1;
  }
  if (caught("fault cages", [&] {
        return seams_frame_cages(oxygen, 3.5, 1, 4, "cutoff", 0, nullptr,
                                 &cageHex, nullptr, nullptr, nullptr) != 0
                   ? 1
                   : 0;
      })) {
    return 1;
  }
  if (caught("fault signature", [&] {
        return seams_frame_signature(oxygen, 3.5, 1, "cutoff", "512",
                                     &faultCages, &faultAtoms) != 0
                   ? 1
                   : 0;
      })) {
    return 1;
  }
  if (caught("fault rdf", [&] {
        return seams_frame_rdf(oxygen, 6.0, 4, 1, 1, faultR, faultG, faultC,
                               nullptr) != 0
                   ? 1
                   : 0;
      })) {
    return 1;
  }
  if (caught("fault rdf2d", [&] {
        return seams_frame_rdf2d(oxygen, 3.5, 0.5, nullptr, faultR, faultG) !=
                       0
                   ? 1
                   : 0;
      })) {
    return 1;
  }
  if (caught("fault prisms", [&] {
        return seams_frame_prisms(oxygen, 3.5, 1, 4, nullptr, faultPrism,
                                  &faultPrismN) != 0
                   ? 1
                   : 0;
      })) {
    return 1;
  }
  if (caught("fault sheet", [&] {
        return seams_frame_monolayer(oxygen, 3.5, 1, 4, 1.0, nullptr,
                                     faultSheet, 8) != 0
                   ? 1
                   : 0;
      })) {
    return 1;
  }
  const long beforeUnknown = fault_live()();
  fault_arm()(1, 2);
  if (seams_frame_chill(oxygen, 3.5, 1, "cutoff", nullptr, nullptr) == 0) {
    return fail("unknown chill");
  }
  fault_arm()(0, 0);
  if (std::string(seams_last_error()) != "unknown error") {
    return fail("unknown text");
  }
  if (fault_live()() != beforeUnknown) {
    return fail("unknown leak");
  }
  }
  if (seams_frame_prisms(oxygen, 3.5, 1, 4, nullptr, nullptr, nullptr) != 0) {
    return fail("prism both null");
  }
  const std::string withSlash = outDir + "/";
  if (seams_frame_rdf2d(oxygen, 3.5, 0.5, withSlash.c_str(), rr, gg) != 0) {
    return fail("rdf2d slash again");
  }
  if (seams_frame_rdf2d(oxygen, 3.5, 0.5, "", rr, gg) != 0) {
    return fail("rdf2d blank dir");
  }
  if (seams_frame_prisms(oxygen, 3.5, 1, 4, withSlash.c_str(), prismWide.data(), &prismCount) != 0) {
    return fail("prism slash again");
  }
  if (seams_frame_prisms(oxygen, 3.5, 1, 4, "", prismWide.data(), &prismCount) != 0) {
    return fail("prism blank dir");
  }
  seams_frame_free(pair);
  {
    const double farXyz[] = {1, 1, 1, 11, 11, 11, 21, 21, 21};
    const double farBox[] = {0, 0, 0, 30, 30, 30, 0, 0, 0};
    seams_frame *far = seams_frame_from_arrays(farXyz, 3, farBox, nullptr, nullptr);
    if (far == nullptr || seams_frame_neighbours(far, 3.5, 1, "mutual") != 0 ||
        seams_frame_neighbour_pairs(far) != 0) {
      seams_frame_free(far);
      return fail("mutual distance bound");
    }
    seams_frame_free(far);
  }
  {
    const double planeXyz[] = {1, 1, 1, 12, 1, 1, 1, 12, 1};
    const double planeBox[] = {0, 0, 0, 30, 30, 30, 0, 0, 0};
    seams_frame *plane =
        seams_frame_from_arrays(planeXyz, 3, planeBox, nullptr, nullptr);
    double planeR[12] = {};
    double planeG[12] = {};
    if (plane == nullptr ||
        seams_frame_rdf2d(plane, 6.0, 0.5, nullptr, planeR, planeG) != 0) {
      seams_frame_free(plane);
      return fail("coplanar rdf2d");
    }
    for (double value : planeG) {
      if (!std::isfinite(value)) {
        seams_frame_free(plane);
        return fail("coplanar rdf2d values");
      }
    }
    seams_frame_free(plane);
  }
  seams_frame_free(blank2);
  seams_frame_free(oxygen);
  std::filesystem::remove_all(out);

  seams_frame_free(molten);
  seams_frame_free(tip);
  seams_frame_free(ice);
  seams_frame_free(guest);
  {
    const double brineXyz[] = {10, 10, 10, 12.8, 10, 10, 9.07, 12.42, 10,
                               9.07, 8.79, 11.62, 9.07, 8.79, 8.38, 12.0, 10, 10};
    const double brineBox[] = {0, 0, 0, 30, 30, 30, 0, 0, 0};
    const int types[] = {1, 1, 1, 1, 1, 3};
    const int mols[] = {1, 2, 3, 4, 5, 6};
    seams_frame *mixed =
        seams_frame_from_species(brineXyz, 6, brineBox, types, nullptr, mols);
    if (mixed == nullptr) {
      return fail("brine frame");
    }
    int mixedLabels[6] = {};
    if (seams_frame_chill_plus(mixed, 3.5, 0, "cutoff", mixedLabels, nullptr) !=
        0) {
      seams_frame_free(mixed);
      return fail("brine untyped chill");
    }
    if (mixedLabels[5] == 6) {
      seams_frame_free(mixed);
      return fail("brine ion was already unclassified");
    }
    if (seams_frame_select(mixed, 1, nullptr, 0) != 0) {
      seams_frame_free(mixed);
      return fail("brine select");
    }
    int selected[6] = {};
    if (seams_frame_chill_plus(mixed, 3.5, 1, "cutoff", selected, nullptr) != 0) {
      seams_frame_free(mixed);
      return fail("brine selected chill");
    }
    if (selected[5] != 6) {
      seams_frame_free(mixed);
      return fail("brine ion stayed labeled");
    }
    int ringCounts[7] = {};
    if (seams_frame_rings(mixed, 3.5, 1, "cutoff", 6, ringCounts, 7) != 0) {
      seams_frame_free(mixed);
      return fail("brine rings");
    }
    int hex = 0;
    int cubicCages = 0;
    int water = 0;
    int both = 0;
    if (seams_frame_cages(mixed, 3.5, 1, 4, "cutoff", 0, nullptr, &hex, &cubicCages,
                          &water, &both) != 0) {
      seams_frame_free(mixed);
      return fail("brine cages");
    }
    if (hex + cubicCages + water + both != 5) {
      seams_frame_free(mixed);
      return fail("brine cages counted the ion");
    }
    if (seams_frame_select(mixed, 9, nullptr, 0) == 0) {
      seams_frame_free(mixed);
      return fail("brine missing type was accepted");
    }
    if (std::string(seams_last_error()) != "atom type 9 is not in the frame") {
      seams_frame_free(mixed);
      return fail("brine missing type text");
    }
    const int absentMol = 99;
    if (seams_frame_select(mixed, 1, &absentMol, 1) == 0) {
      seams_frame_free(mixed);
      return fail("brine missing molecule was accepted");
    }
    seams_frame_free(mixed);
  }
  {
    seams_frame *fresh = seams_frame_read(tip4p.c_str(), 1, 2);
    if (fresh == nullptr) {
      return fail("reread tip4p");
    }
    const int nFresh = seams_frame_nop(fresh);
    std::vector<double> moved(static_cast<std::size_t>(nFresh) * 3);
    const double *src = seams_frame_positions(fresh);
    for (int i = 0; i < nFresh * 3; ++i) {
      moved[static_cast<std::size_t>(i)] = src[i];
    }
    moved[0] += 1.0e-4;
    double freshBox[9];
    seams_frame_box(fresh, freshBox);
    if (seams_frame_set_positions(fresh, moved.data(), nFresh, freshBox) != 0) {
      seams_frame_free(fresh);
      return fail("set positions");
    }
    if (seams_frame_neighbours(fresh, 3.5, 2, "cutoff") != 0) {
      seams_frame_free(fresh);
      return fail("neighbours");
    }
    if (seams_frame_neighbour_pairs(fresh) <= 0) {
      seams_frame_free(fresh);
      return fail("neighbour pairs");
    }
    std::vector<double> freshQl(static_cast<std::size_t>(nFresh));
    std::vector<double> bar(static_cast<std::size_t>(nFresh));
    if (seams_frame_steinhardt(fresh, 3.5, 2, 6, freshQl.data(), bar.data()) != 0) {
      seams_frame_free(fresh);
      return fail("steinhardt");
    }
    if (seams_frame_steinhardt_graph(fresh, 3.5, 2, 6, "cutoff", freshQl.data(), bar.data()) != 0) {
      seams_frame_free(fresh);
      return fail("steinhardt graph");
    }
    std::vector<double> wl(static_cast<std::size_t>(nFresh));
    if (seams_frame_wl(fresh, 3.5, 2, 4, wl.data()) != 0) {
      seams_frame_free(fresh);
      return fail("wl");
    }
    std::vector<double> q4(static_cast<std::size_t>(nFresh));
    std::vector<double> q6(static_cast<std::size_t>(nFresh));
    if (seams_frame_q4_q6(fresh, 3.5, 2, q4.data(), q6.data(), nullptr, nullptr) != 0) {
      seams_frame_free(fresh);
      return fail("q4 q6");
    }
    std::vector<int> iceLabels(static_cast<std::size_t>(nFresh));
    int fused[9] = {};
    if (seams_frame_q4_q6_chill(fresh, 3.5, 2, q4.data(), q6.data(), nullptr, nullptr,
                               iceLabels.data(), fused) != 0) {
      seams_frame_free(fresh);
      return fail("q4 q6 chill");
    }
    std::vector<int> ringCounts(9);
    seams_phase_times times{};
    if (seams_frame_chill_rings(fresh, 3.5, 2, "cutoff", 6, 0.0, iceLabels.data(), fused,
                               ringCounts.data(), 9, &times) != 0) {
      seams_frame_free(fresh);
      return fail("chill rings");
    }
    std::vector<int> shell(static_cast<std::size_t>(nFresh));
    if (seams_frame_symmetry(fresh, 0.12, shell.data()) != 0) {
      seams_frame_free(fresh);
      return fail("symmetry");
    }
    double tau = 0.0;
    if (seams_frames_jump_tau90(fresh, fresh, 1.0, 2, 1, &tau) != 0) {
      seams_frame_free(fresh);
      return fail("jump tau");
    }
    seams_frame_free(fresh);
  }
  std::cout << "seams_lang ok " << seams_lang_version() << "\n";
  const std::string shPath = root + "/input/traj/genice_sH.lammpstrj";
  seams_frame *sh = seams_frame_read(shPath.c_str(), 1, 1);
  if (sh == nullptr || seams_frame_nop(sh) != 136) {
    return fail("sH frame");
  }
  double shBox[9] = {};
  seams_frame_box(sh, shBox);
  seams_frame_free(sh);
  if (std::abs(shBox[3] - 24.8424228) > 1e-4 || std::abs(shBox[4] - 10.7570846) > 1e-4 ||
      std::abs(shBox[5] - 20.1312910) > 1e-4 || std::abs(shBox[6] + 6.2106057) > 1e-4 ||
      std::abs(shBox[7]) > 1e-6 || std::abs(shBox[8]) > 1e-6) {
    std::cerr << "sH box";
    for (double value : shBox) {
      std::cerr << " " << value;
    }
    std::cerr << "\n";
    return fail("sH cell");
  }
  const double heavy[] = {0, 0, 0, 2.8, 0, 0};
  const double hxyz[] = {1, 0, 0};
  const int heavyTypes[] = {8, 8};
  const int heavyMols[] = {1, 2};
  const int hydrogenTypes[] = {2};
  const int hydrogenMols[] = {1};
  const double linBox[] = {0, 0, 0, 20, 20, 20, 0, 0, 0};
  seams_frame *hframe = seams_frame_from_hydrogen_cloud(
      heavy, 2, heavyTypes, nullptr, heavyMols, hxyz, 1, hydrogenTypes, nullptr, hydrogenMols,
      linBox);
  if (hframe == nullptr || seams_frame_nop(hframe) != 3) {
    return fail("hydrogen cloud");
  }
  const int donors[] = {8};
  const int acceptors[] = {10};
  int n_dir = -1;
  if (seams_frame_directed_hbonds(hframe, 3.5, 2.42, 30.0, donors, 1, acceptors, 0, 2, 1, nullptr,
                                  nullptr, nullptr, nullptr, nullptr, nullptr, 0, &n_dir) != 0 ||
      n_dir != 0) {
    seams_frame_free(hframe);
    return fail("empty acceptor set");
  }
  uint32_t donorIdx = 0;
  uint32_t hyd = 0;
  uint32_t acc = 0;
  double dah = 0;
  double dda = 0;
  double ang = 0;
  n_dir = -1;
  if (seams_frame_directed_hbonds(hframe, 3.5, 2.42, 30.0, donors, 1, donors, 1, 2, 1, &donorIdx,
                                  &hyd, &acc, &dah, &dda, &ang, 1, &n_dir) != 0 ||
      n_dir != 1 || donorIdx != 0 || hyd != 2 || acc != 1 || std::abs(dda - 2.8) > 1e-12 ||
      std::abs(dah - 1.8) > 1e-12 || std::abs(ang) > 1e-8) {
    seams_frame_free(hframe);
    return fail("directed linear bond");
  }
  seams_frame_free(hframe);
  if (seams_frame_from_hydrogen_cloud(nullptr, 2, heavyTypes, nullptr, heavyMols, hxyz, 1,
                                      hydrogenTypes, nullptr, hydrogenMols, linBox) != nullptr) {
    return fail("null heavy cloud");
  }
  return 0;
}
