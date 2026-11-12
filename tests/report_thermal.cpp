//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// Reads frames written by scripts/mw_frames.py. Prints ice fractions, the
// stacking sequence, the eclipsed-cutoff sweep, and a per-atom file that
// the PTM pass joins to.

#include <local_symmetry.hpp>
#include <seams_c_api.h>

#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Frame {
  std::string name;
  double temperature = 0;
  std::string role;
  Eigen::Matrix3d box = Eigen::Matrix3d::Zero();
  Eigen::MatrixXd xyz;
  std::vector<int> bilayer;
  std::vector<int> core;
};

const char *iceName(int label) {
  if (label == 0) {
    return "Ih";
  }
  if (label == 1) {
    return "Ic";
  }
  if (label == 2) {
    return "Isd";
  }
  return "other";
}

const char *chillName(int label) {
  const char *names[] = {"cubic", "hexagonal", "water", "interfacial", "clathrate",
                         "interClathrate", "unclassified", "reCubic", "reHex"};
  if (label >= 0 && label <= 8) {
    return names[label];
  }
  return "bad";
}

int chillFamily(int label) {
  if (label == 0 || label == 7) {
    return 1;
  }
  if (label == 1 || label == 8) {
    return 0;
  }
  return -1;
}

bool loadFrame(const std::string &path, Frame &frame) {
  std::ifstream in(path);
  if (!in) {
    return false;
  }
  std::string tag;
  int n = 0;
  in >> tag >> frame.name;
  in >> tag >> frame.temperature;
  in >> tag >> frame.role;
  in >> tag >> n;
  in >> tag;
  for (int r = 0; r < 3; r++) {
    for (int c = 0; c < 3; c++) {
      in >> frame.box(r, c);
    }
  }
  in >> tag;
  frame.xyz.resize(n, 3);
  frame.bilayer.assign(static_cast<size_t>(n), -1);
  frame.core.assign(static_cast<size_t>(n), 0);
  for (int i = 0; i < n; i++) {
    in >> frame.xyz(i, 0) >> frame.xyz(i, 1) >> frame.xyz(i, 2) >> frame.bilayer[static_cast<size_t>(i)] >>
        frame.core[static_cast<size_t>(i)];
  }
  return static_cast<bool>(in);
}

void histogram(const std::vector<int> &labels, const char *(*name)(int)) {
  std::map<std::string, int> hist;
  for (int label : labels) {
    hist[name(label)] += 1;
  }
  for (const auto &item : hist) {
    std::cout << " " << item.first << "=" << item.second;
  }
}

int geometryOf(const Eigen::MatrixXd &shell, double &rIh, double &rIc, double &rIsd) {
  rIh = symm::iceTemplateRmsd(shell, symm::Kind::iceIh);
  rIc = symm::iceTemplateRmsd(shell, symm::Kind::iceIc);
  rIsd = symm::iceTemplateRmsd(shell, symm::Kind::iceIsd);
  const double r[3] = {rIh, rIc, rIsd};
  int best = -1;
  double bestR = 1e300;
  double second = 1e300;
  for (int k = 0; k < 3; k++) {
    if (!(r[k] >= 0.0)) {
      continue;
    }
    if (r[k] < bestR) {
      second = bestR;
      bestR = r[k];
      best = k;
    } else if (r[k] < second) {
      second = r[k];
    }
  }
  if (best < 0 || !(bestR < 0.55)) {
    return -1;
  }
  if (second - bestR < 0.03) {
    return -2;
  }
  return best;
}

void analyze(const Frame &frame, const std::string &path) {
  const int n = static_cast<int>(frame.xyz.rows());
  const double cutoffs[] = {0.75, 0.85, 0.90, 0.95, 1.05};
  std::vector<std::vector<int>> byCut(5);
  for (int k = 0; k < 5; k++) {
    symm::classifyIce(frame.xyz, frame.box, byCut[static_cast<size_t>(k)], cutoffs[k]);
  }
  const std::vector<int> &kernel = byCut[2];
  std::vector<double> raw(static_cast<size_t>(3 * n));
  for (int i = 0; i < n; i++) {
    raw[static_cast<size_t>(3 * i)] = frame.xyz(i, 0);
    raw[static_cast<size_t>(3 * i + 1)] = frame.xyz(i, 1);
    raw[static_cast<size_t>(3 * i + 2)] = frame.xyz(i, 2);
  }
  const double box9[9] = {0,
                          0,
                          0,
                          frame.box(0, 0),
                          frame.box(1, 1),
                          frame.box(2, 2),
                          frame.box(1, 0),
                          frame.box(2, 0),
                          frame.box(2, 1)};
  std::vector<int> chill(static_cast<size_t>(n), -1);
  if (seams_chill_plus_graph_lammps(raw.data(), n, box9, 3.5, "cutoff", chill.data()) != 0) {
    std::cout << frame.name << " chill failed\n";
    return;
  }

  std::cout << "frame " << frame.name << " T=" << frame.temperature << " role=" << frame.role
            << " n=" << n << "\n";
  std::cout << "kernel";
  histogram(kernel, iceName);
  std::cout << "\nchill";
  histogram(chill, chillName);
  std::cout << "\n";

  std::map<int, std::array<int, 4>> layers;
  int coreIce = 0;
  int coreN = 0;
  int outerIce = 0;
  int outerN = 0;
  int agree = 0;
  int kernelRight = 0;
  int chillRight = 0;
  int bothRight = 0;
  int neither = 0;
  int ambiguous = 0;
  int contested = 0;
  std::ofstream labels(path + ".labels");
  labels << "index kernel chill geom rmsd_ih rmsd_ic rmsd_isd bilayer core "
            "c0 c1 c2 c3\n";
  for (int i = 0; i < n; i++) {
    const Eigen::MatrixXd shell = symm::iceShell(frame.xyz, frame.box, i);
    double rIh = -1;
    double rIc = -1;
    double rIsd = -1;
    double cosines[4] = {-2, -2, -2, -2};
    const int geom = shell.rows() == 17 ? geometryOf(shell, rIh, rIc, rIsd) : -1;
    if (shell.rows() == 17) {
      (void)symm::bondCosines(shell, cosines);
    }
    const int kLab = kernel[static_cast<size_t>(i)];
    const int cLab = chillFamily(chill[static_cast<size_t>(i)]);
    labels << i << " " << kLab << " " << chill[static_cast<size_t>(i)] << " " << geom << " " << rIh
           << " " << rIc << " " << rIsd << " " << frame.bilayer[static_cast<size_t>(i)] << " "
           << frame.core[static_cast<size_t>(i)];
    for (double c : cosines) {
      labels << " " << c;
    }
    labels << "\n";
    if (kLab == cLab && kLab >= 0) {
      agree++;
    }
    const bool clash = kLab != cLab;
    if (geom == -2) {
      ambiguous++;
    } else if (clash) {
      contested++;
      const bool kOk = kLab == geom;
      const bool cOk = cLab == geom;
      if (kOk && cOk) {
        bothRight++;
      } else if (kOk) {
        kernelRight++;
      } else if (cOk) {
        chillRight++;
      } else {
        neither++;
      }
    }
    const int layer = frame.bilayer[static_cast<size_t>(i)];
    if (layer >= 0) {
      int slot = 3;
      if (kLab >= 0 && kLab <= 2) {
        slot = kLab;
      }
      layers[layer][static_cast<size_t>(slot)] += 1;
    }
    if (frame.core[static_cast<size_t>(i)]) {
      coreN++;
      if (kLab >= 0) {
        coreIce++;
      }
    } else if (frame.role == "nucleus") {
      outerN++;
      if (kLab >= 0) {
        outerIce++;
      }
    }
  }
  if (!layers.empty()) {
    std::cout << "sequence";
    for (const auto &item : layers) {
      const auto &c = item.second;
      const int labelled = c[0] + c[1] + c[2];
      int winner = 3;
      int best = c[3];
      for (int s = 0; s < 3; s++) {
        if (c[static_cast<size_t>(s)] > best) {
          best = c[static_cast<size_t>(s)];
          winner = s;
        }
      }
      const char *token = "x";
      if (winner == 0) {
        token = "h";
      } else if (winner == 1) {
        token = "c";
      } else if (winner == 2) {
        token = "d";
      }
      const bool mixed = labelled > 0 && best * 4 < labelled * 3;
      std::cout << " " << item.first << ":" << (mixed ? "m" : "") << token << "(Ih=" << c[0]
                << ",Ic=" << c[1] << ",Isd=" << c[2] << ",other=" << c[3] << ")";
    }
    std::cout << "\n";
  }
  if (frame.role == "nucleus") {
    std::cout << "core_ice " << coreIce << " of " << coreN << " outer_ice " << outerIce << " of "
              << outerN << "\n";
  }
  std::cout << "kernel_chill_same_family " << agree << " contested " << contested
            << " kernel_matches_shell " << kernelRight << " chill_matches_shell " << chillRight
            << " both " << bothRight << " neither " << neither << " ambiguous_shell " << ambiguous
            << "\n";
  for (int k = 0; k < 5; k++) {
    int changed = 0;
    for (int i = 0; i < n; i++) {
      if (byCut[static_cast<size_t>(k)][static_cast<size_t>(i)] != kernel[static_cast<size_t>(i)]) {
        changed++;
      }
    }
    std::cout << "cutoff " << cutoffs[k] << " changed " << changed << " of " << n << "\n";
  }
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "usage: report_thermal frame [frame...]\n";
    return 2;
  }
  for (int i = 1; i < argc; i++) {
    Frame frame;
    if (!loadFrame(argv[i], frame)) {
      std::cerr << "unreadable " << argv[i] << "\n";
      return 1;
    }
    analyze(frame, argv[i]);
  }
  return 0;
}
