#include "seams_lang.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

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
  if (count_of(counts, 5) != 12 || count_of(counts, 2) != 238 ||
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
  if (std::abs(mean - 0.389362) > 1e-4) {
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
  for (int i = 3; i <= 8; ++i) {
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
    FILE *pipe = popen(cmd.c_str(), "r");
    if (pipe == nullptr) {
      return fail("cli pipe");
    }
    std::string line;
    char buf[512];
    while (fgets(buf, sizeof buf, pipe) != nullptr) {
      line += buf;
    }
    const int rc = pclose(pipe);
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

  seams_frame_free(molten);
  seams_frame_free(tip);
  seams_frame_free(ice);
  seams_frame_free(guest);
  std::cout << "seams_lang ok " << seams_lang_version() << "\n";
  return 0;
}
