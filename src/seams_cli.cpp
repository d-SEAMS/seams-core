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

void printCounts(const Cloud &cloud) {
  std::map<std::string, int> hist;
  for (const auto &pt : cloud.pts) {
    hist[iceName(pt.iceType)]++;
  }
  std::cout << "nop " << cloud.nop;
  for (const auto &[name, n] : hist) {
    if (n > 0) {
      std::cout << " " << name << " " << n;
    }
  }
  std::cout << "\n";
}

int cmdRead(Cloud &cloud) {
  const auto box = cloud.box;
  std::cout << "nop " << cloud.nop << " frame " << cloud.currentFrame
            << " box " << box[0] << " " << box[1] << " " << box[2] << "\n";
  return 0;
}

int cmdChillPlus(Cloud &cloud, double cutoff, int typeI) {
  const int typ = typeOf(cloud, typeI);
  auto nList = nneigh::neighListO(cutoff, cloud, typ);
  chill::getCorrelPlus(cloud, nList, false);
  chill::getIceTypePlusNoPrint(cloud, nList, false);
  printCounts(cloud);
  return 0;
}

int cmdChill(Cloud &cloud, double cutoff, int typeI) {
  const int typ = typeOf(cloud, typeI);
  auto nList = nneigh::neighListO(cutoff, cloud, typ);
  chill::getCorrel(cloud, nList, false);
  chill::getIceTypeNoPrint(cloud, nList, false);
  printCounts(cloud);
  return 0;
}

int cmdCages(Cloud &cloud, double cutoff, int typeI, int k) {
  const int typ = typeOf(cloud, typeI);
  const double cand = cutoff + 1.5;
  auto mutual = nneigh::kNearestNeighbourList(cloud, k, cand, typ, true);
  auto uni = nneigh::kNearestNeighbourList(cloud, k, cand, typ, false);
  auto idxS = nneigh::neighbourListByIndex(cloud, mutual);
  auto idxU = nneigh::neighbourListByIndex(cloud, uni);
  auto ringsS = primitive::ringNetwork(idxS, 6);
  auto ringsU = primitive::ringNetwork(idxU, 6);
  std::vector<std::vector<int>> sixS;
  std::vector<std::vector<int>> sixU;
  for (const auto &r : ringsS) {
    if (r.size() == 6) {
      sixS.push_back(r);
    }
  }
  for (const auto &r : ringsU) {
    if (r.size() == 6) {
      sixU.push_back(r);
    }
  }
  const auto aff = ring::seededCageAffiliation(sixS, idxS, sixU, idxU);
  int ih = 0;
  int ic = 0;
  int water = 0;
  for (size_t i = 0; i < aff.hc.size(); ++i) {
    if (aff.hc[i]) {
      ++ih;
    } else if (aff.ddc[i]) {
      ++ic;
    } else {
      ++water;
    }
  }
  std::cout << "nop " << cloud.nop << " ih " << ih << " ic " << ic << " water "
            << water << "\n";
  return 0;
}

} // namespace

int main(int argc, char *argv[]) {
  cxxopts::Options opt(
      argv[0], "d-SEAMS engine CLI. Lua is the luadseams library; Python is pydseams.");
  opt.add_options()("h,help", "Print help")("v,version", "Print version")(
      "f,frame", "Frame number (1-based)",
      cxxopts::value<int>()->default_value("1"))(
      "t,type", "Atom type (0 guesses oxygen then type 1)",
      cxxopts::value<int>()->default_value("0"))(
      "c,cutoff", "Neighbour cutoff (Angstrom)",
      cxxopts::value<double>()->default_value("3.5"))(
      "k", "k for seeded cages", cxxopts::value<int>()->default_value("4"))(
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
  const int typeI = args["type"].as<int>();
  const double cutoff = args["cutoff"].as<double>();
  const int k = args["k"].as<int>();

  Cloud cloud = load(file, frame, typeI);
  if (cmd == "read") {
    return cmdRead(cloud);
  }
  if (cmd == "chill-plus" || cmd == "chill_plus") {
    return cmdChillPlus(cloud, cutoff, typeI);
  }
  if (cmd == "chill") {
    return cmdChill(cloud, cutoff, typeI);
  }
  if (cmd == "cages") {
    return cmdCages(cloud, cutoff, typeI, k);
  }
  std::cerr << "unknown command: " << cmd << "\n";
  return 2;
}
