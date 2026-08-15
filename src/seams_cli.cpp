//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <bop.hpp>
#include <cage_affiliation.hpp>
#include <cxxopts.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

const char *iceName(molSys::atom_state_type s) {
  switch (s) {
  case molSys::atom_state_type::cubic:
    return "cubic";
  case molSys::atom_state_type::hexagonal:
    return "hexagonal";
  case molSys::atom_state_type::water:
    return "water";
  case molSys::atom_state_type::interfacial:
    return "interfacial";
  case molSys::atom_state_type::clathrate:
    return "clathrate";
  case molSys::atom_state_type::interClathrate:
    return "interClathrate";
  case molSys::atom_state_type::reCubic:
    return "reCubic";
  case molSys::atom_state_type::reHex:
    return "reHex";
  case molSys::atom_state_type::unclassified:
    break;
  }
  return "unclassified";
}

Cloud load(const std::string &path, int frame, int typeI) {
  Cloud cloud;
  const auto ext = path.substr(path.find_last_of('.') + 1);
  if (ext == "xyz") {
    return sinp::readXYZ(path);
  }
#ifdef SEAMS_HAS_READCON
  if (ext == "con") {
    return sinp::readCon(path, frame, cloud);
  }
#endif
#ifdef SEAMS_HAS_CHEMFILES
  if (ext == "pdb" || ext == "gro" || ext == "dcd") {
    return sinp::readChemfiles(path, frame, cloud, typeI > 0 ? typeI : -1);
  }
#endif
  if (typeI > 0) {
    return sinp::readLammpsTrjO(path, frame, cloud, typeI);
  }
  cloud = sinp::readLammpsTrjO(path, frame, cloud, 2);
  if (cloud.nop > 0) {
    return cloud;
  }
  Cloud again;
  return sinp::readLammpsTrjO(path, frame, again, 1);
}

int typeOf(const Cloud &cloud, int requested) {
  if (requested > 0) {
    return requested;
  }
  if (cloud.nop == 0) {
    return 1;
  }
  return cloud.pts[0].type;
}

void printCounts(std::ostream &os, const Cloud &cloud) {
  std::map<std::string, int> hist;
  for (const auto &pt : cloud.pts) {
    hist[iceName(pt.iceType)]++;
  }
  os << "nop " << cloud.nop;
  for (const auto &[name, n] : hist) {
    if (n > 0) {
      os << " " << name << " " << n;
    }
  }
  os << "\n";
}

int cmdRead(std::ostream &os, Cloud &cloud) {
  const auto box = cloud.box;
  os << "nop " << cloud.nop << " frame " << cloud.currentFrame << " box "
     << box[0] << " " << box[1] << " " << box[2] << "\n";
  return 0;
}

int cmdChillPlus(std::ostream &os, Cloud &cloud, double cutoff, int typeI) {
  const int typ = typeOf(cloud, typeI);
  auto nList = nneigh::neighListO(cutoff, cloud, typ);
  chill::getCorrelPlus(cloud, nList, false);
  chill::getIceTypePlusNoPrint(cloud, nList, false);
  printCounts(os, cloud);
  return 0;
}

int cmdChill(std::ostream &os, Cloud &cloud, double cutoff, int typeI) {
  const int typ = typeOf(cloud, typeI);
  auto nList = nneigh::neighListO(cutoff, cloud, typ);
  chill::getCorrel(cloud, nList, false);
  chill::getIceTypeNoPrint(cloud, nList, false);
  printCounts(os, cloud);
  return 0;
}

int cmdCages(std::ostream &os, Cloud &cloud, double cutoff, int typeI, int k,
             const std::string &graphName) {
  const int typ = typeOf(cloud, typeI);
  const double cand = cutoff + 1.5;

  auto sixOf = [](const std::vector<std::vector<int>> &rings) {
    std::vector<std::vector<int>> six;
    for (const auto &r : rings) {
      if (r.size() == 6) {
        six.push_back(r);
      }
    }
    return six;
  };

  int ih = 0;
  int ic = 0;
  int water = 0;
  const auto tallyAtoms = [&](const std::vector<bool> &hc,
                              const std::vector<bool> &ddc) {
    ih = ic = water = 0;
    const int n = static_cast<int>(hc.size());
    for (int i = 0; i < n; ++i) {
      if (hc[static_cast<std::size_t>(i)]) {
        ++ih;
      } else if (ddc[static_cast<std::size_t>(i)]) {
        ++ic;
      } else {
        ++water;
      }
    }
  };

  if (graphName == "seeded") {
    auto mutual = nneigh::kNearestNeighbourList(cloud, k, cand, typ, true);
    auto uni = nneigh::kNearestNeighbourList(cloud, k, cand, typ, false);
    auto idxS = nneigh::neighbourListByIndex(cloud, mutual);
    auto idxU = nneigh::neighbourListByIndex(cloud, uni);
    auto sixS = sixOf(primitive::ringNetwork(idxS, 6));
    auto sixU = sixOf(primitive::ringNetwork(idxU, 6));
    const auto aff = ring::seededCageAffiliation(sixS, idxS, sixU, idxU);
    tallyAtoms(aff.hc, aff.ddc);
  } else {
    const auto graph = nneigh::bondGraphFromName(graphName);
    std::vector<std::vector<int>> nList;
    if (graph == nneigh::BondGraph::Cutoff) {
      nList = nneigh::neighListO(cutoff, cloud, typ);
    } else {
      const bool mutual = graph == nneigh::BondGraph::KnnMutual;
      nList = nneigh::kNearestNeighbourList(cloud, k, cand, typ, mutual);
    }
    auto idx = nneigh::neighbourListByIndex(cloud, nList);
    auto six = sixOf(primitive::ringNetwork(idx, 6));
    const auto aff = ring::cageAffiliation(six, idx);
    // cageAffiliation is per-ring; map to atoms
    std::vector<bool> hc(static_cast<std::size_t>(cloud.nop), false);
    std::vector<bool> ddc(static_cast<std::size_t>(cloud.nop), false);
    for (std::size_t r = 0; r < six.size(); ++r) {
      for (const int a : six[r]) {
        if (a >= 0 && a < cloud.nop) {
          hc[static_cast<std::size_t>(a)] =
              hc[static_cast<std::size_t>(a)] || aff.hc[r];
          ddc[static_cast<std::size_t>(a)] =
              ddc[static_cast<std::size_t>(a)] || aff.ddc[r];
        }
      }
    }
    tallyAtoms(hc, ddc);
  }
  os << "nop " << cloud.nop << " graph " << graphName << " ih " << ih
     << " ic " << ic << " water " << water << "\n";
  return 0;
}

} // namespace

int main(int argc, char *argv[]) {
  cxxopts::Options opt(
      argv[0], "d-SEAMS engine CLI. Lua is the luadseams library; Python is pydseams.");
  opt.add_options()("h,help", "Print help")("v,version", "Print version")(
      "f,frame", "First frame (1-based)",
      cxxopts::value<int>()->default_value("1"))(
      "last", "Last frame (inclusive). Omit for a single --frame.",
      cxxopts::value<int>()->default_value("0"))(
      "j,jobs", "Parallel frame workers (OpenMP). 1 is serial.",
      cxxopts::value<int>()->default_value("1"))(
      "t,type", "Atom type (0 guesses oxygen then type 1)",
      cxxopts::value<int>()->default_value("0"))(
      "c,cutoff", "Neighbour cutoff (Angstrom)",
      cxxopts::value<double>()->default_value("3.5"))(
      "k", "k for knn / seeded cages", cxxopts::value<int>()->default_value("4"))(
      "graph",
      "Bond graph for cages: cutoff | knn | knn-union | seeded",
      cxxopts::value<std::string>()->default_value("seeded"))(
      "command", "read | chill | chill-plus | cages",
      cxxopts::value<std::string>())("file", "Trajectory file",
                                     cxxopts::value<std::string>());
  opt.parse_positional({"command", "file"});
  opt.positional_help("COMMAND FILE");

  cxxopts::ParseResult args;
  try {
    args = opt.parse(argc, argv);
  } catch (const cxxopts::OptionException &e) {
    std::cerr << e.what() << "\n";
    return 2;
  }
  if (args.count("help")) {
    std::cout << opt.help() << "\n";
    return 0;
  }
  if (args.count("version")) {
#ifdef SEAMS_VERSION
    std::cout << "seams " << SEAMS_VERSION << "\n";
#else
    std::cout << "seams\n";
#endif
    return 0;
  }
  if (!args.count("command") || !args.count("file")) {
    std::cerr << opt.help() << "\n";
    return 2;
  }
  const std::string cmd = args["command"].as<std::string>();
  const std::string file = args["file"].as<std::string>();
  const int frame = args["frame"].as<int>();
  const int last = args["last"].as<int>();
  const int jobs = args["jobs"].as<int>();
  const int typeI = args["type"].as<int>();
  const double cutoff = args["cutoff"].as<double>();
  const int k = args["k"].as<int>();
  const std::string graph = args["graph"].as<std::string>();

  auto runOne = [&](std::ostream &os, Cloud &cloud) {
    if (cmd == "read") {
      return cmdRead(os, cloud);
    }
    if (cmd == "chill-plus" || cmd == "chill_plus") {
      return cmdChillPlus(os, cloud, cutoff, typeI);
    }
    if (cmd == "chill") {
      return cmdChill(os, cloud, cutoff, typeI);
    }
    if (cmd == "cages") {
      try {
        return cmdCages(os, cloud, cutoff, typeI, k, graph);
      } catch (const std::exception &e) {
        os << e.what() << "\n";
        return 2;
      }
    }
    os << "unknown command: " << cmd << "\n";
    return 2;
  };

  if (last <= 0 || last == frame) {
    Cloud cloud = load(file, frame, typeI);
    return runOne(std::cout, cloud);
  }

  std::mutex outMu;
  int rc = 0;
  const int typeFilter = typeI > 0 ? typeI : 0;
  sinp::forEachLammpsFrame(
      file, frame, last, typeFilter,
      [&](int /*fr*/, Cloud &cloud) {
        if (typeFilter <= 0 && cloud.nop == 0) {
          cloud = load(file, cloud.currentFrame, typeI);
        }
        std::ostringstream line;
        const int one = runOne(line, cloud);
        {
          std::lock_guard<std::mutex> lock(outMu);
          std::cout << line.str();
          if (one != 0) {
            rc = one;
          }
        }
      },
      jobs);
  return rc;
}
