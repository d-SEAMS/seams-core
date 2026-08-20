//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <argum.h>
#include <bond.hpp>
#include <bop.hpp>
#include <cage_affiliation.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <generic.hpp>
#include <neighbours.hpp>
#include <cluster.hpp>
#include <density.hpp>
#include <rdf.hpp>
#include <seams_config.hpp>
#include <seams_input.hpp>
#include <site.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using namespace Argum;

Colorizer colorizer;

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
    const int filter = typeI < 0 ? -1 : (typeI > 0 ? typeI : -1);
    return sinp::readChemfiles(path, frame, cloud, filter);
  }
#endif
  if (typeI < 0) {
    return sinp::readLammpsTrj(path, frame, cloud);
  }
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

std::string_view trimView(std::string_view s) {
  const auto first = s.find_first_not_of(" \t");
  if (first == std::string_view::npos) {
    return {};
  }
  const auto last = s.find_last_not_of(" \t");
  return s.substr(first, last - first + 1);
}

template <typename Fn>
bool forEachInputFrame(const std::string &path, int first, int last,
                       int typeFilter, int nThreads, Fn &&fn) {
  const auto dot = path.find_last_of('.');
  const std::string ext = dot == std::string::npos ? "" : path.substr(dot + 1);
  if (ext == "lammpstrj" || ext == "dump" || ext == "trj") {
    sinp::forEachLammpsFrame(path, first, last, typeFilter,
                             std::forward<Fn>(fn), nThreads);
    return true;
  }

  if (ext == "xyz") {
    if (first <= 1 && (last <= 0 || last >= 1)) {
      Cloud cloud = load(path, 1, typeFilter);
      fn(1, cloud);
    }
    return true;
  }

  if (ext == "con" || ext == "pdb" || ext == "gro" || ext == "dcd") {
    const int stop = last <= 0 ? first : last;
    for (int frame = std::max(1, first); frame <= stop; ++frame) {
      Cloud cloud = load(path, frame, typeFilter);
      if (cloud.nop == 0 && frame != std::max(1, first)) {
        break;
      }
      fn(frame, cloud);
      if (cloud.nop == 0) {
        break;
      }
    }
    return true;
  }

  sinp::forEachLammpsFrame(path, first, last, typeFilter,
                           std::forward<Fn>(fn), nThreads);
  return true;
}

std::string jsonEscape(std::string_view value) {
  std::string escaped;
  escaped.reserve(value.size());
  for (const unsigned char ch : value) {
    switch (ch) {
    case '\\':
      escaped += "\\\\";
      break;
    case '"':
      escaped += "\\\"";
      break;
    case '\n':
      escaped += "\\n";
      break;
    case '\r':
      escaped += "\\r";
      break;
    case '\t':
      escaped += "\\t";
      break;
    default:
      if (ch < 0x20) {
        constexpr char hex[] = "0123456789abcdef";
        escaped += "\\u00";
        escaped += hex[(ch >> 4) & 0xf];
        escaped += hex[ch & 0xf];
      } else {
        escaped.push_back(static_cast<char>(ch));
      }
    }
  }
  return escaped;
}

void emitResult(std::ostream &os, std::string_view format,
                std::string_view command, int frame, int status,
                std::string_view text) {
  if (format != "json") {
    os << text;
    return;
  }
  os << "{\"schema\":\"dseams.cli/v1\",\"command\":\""
     << jsonEscape(command) << "\",\"frame\":" << frame
     << ",\"status\":" << status << ",\"text\":\""
     << jsonEscape(text) << "\"}\n";
}

bool parseTypePair(std::string_view value, int &typeI, int &typeJ) {
  const auto comma = value.find(',');
  if (comma == std::string_view::npos) {
    return false;
  }
  if (value.find(',', comma + 1) != std::string_view::npos) {
    return false;
  }
  const auto left = trimView(value.substr(0, comma));
  const auto right = trimView(value.substr(comma + 1));
  if (left.empty() || right.empty()) {
    return false;
  }
  try {
    typeI = parseIntegral<int>(left);
    typeJ = parseIntegral<int>(right);
  } catch (const ParsingException &) {
    return false;
  }
  return true;
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

std::string iceColor(std::string_view name) {
  if (name == "cubic" || name == "reCubic") {
    return colorizer.longOption(name);
  }
  if (name == "hexagonal" || name == "reHex") {
    return colorizer.shortOption(name);
  }
  if (name == "interfacial") {
    return colorizer.warning(name);
  }
  if (name == "clathrate" || name == "interClathrate") {
    return colorizer.progName(name);
  }
  if (name == "water") {
    return std::string(name);
  }
  return colorizer.error(name);
}

void printCounts(std::ostream &os, const Cloud &cloud) {
  std::map<std::string, int> hist;
  for (const auto &pt : cloud.pts) {
    hist[iceName(pt.iceType)]++;
  }
  os << colorizer.heading("nop") << " " << cloud.nop;
  for (const auto &[name, n] : hist) {
    if (n > 0) {
      os << " " << iceColor(name) << " " << n;
    }
  }
  os << "\n";
}

void printFeatures(std::ostream &os) {
  os << colorizer.heading("Compile-time features:") << "\n";
  const auto line = [&](const char *name, bool on) {
    os << name << ": "
       << (on ? colorizer.warning("enabled") : colorizer.error("disabled"))
       << "\n";
  };
#ifdef SEAMS_HAS_OPENMP
  line("OpenMP", true);
#else
  line("OpenMP", false);
#endif
#ifdef SEAMS_HAS_MPI
  line("MPI", true);
#else
  line("MPI", false);
#endif
#ifdef SEAMS_HAS_HWY
  line("Highway SIMD", true);
#else
  line("Highway SIMD", false);
#endif
#ifdef SEAMS_HAS_VESIN
  line("vesin neighbours", true);
#else
  line("vesin neighbours", false);
#endif
#ifdef SEAMS_HAS_LINKCELL
  line("linkcell k-nearest", true);
#else
  line("linkcell k-nearest", false);
#endif
#ifdef SEAMS_HAS_CHEMFILES
  line("chemfiles", true);
#else
  line("chemfiles", false);
#endif
#ifdef SEAMS_HAS_READCON
  line("readcon-core", true);
#else
  line("readcon-core", false);
#endif
#ifdef SEAMS_HAS_IRA
  line("IRA/SOFI", true);
#else
  line("IRA/SOFI", false);
#endif
#ifdef SEAMS_HAS_SPHERICART
  line("sphericart", true);
#else
  line("sphericart", false);
#endif
#ifdef SEAMS_HAS_NAUTY
  line("nauty", true);
#else
  line("nauty", false);
#endif
}

int cmdRead(std::ostream &os, Cloud &cloud) {
  os << colorizer.heading("nop") << " " << cloud.nop << " "
     << colorizer.longOption("frame") << " " << cloud.currentFrame << " "
     << colorizer.shortOption("box") << " " << gen::formatDumpBox(cloud.box)
     << "\n";
  return 0;
}

site::Kind kindFromName(std::string_view name) {
  if (name == "unspecified") {
    return site::Kind::unspecified;
  }
  if (name == "cationHead") {
    return site::Kind::cationHead;
  }
  if (name == "anion") {
    return site::Kind::anion;
  }
  if (name == "tail") {
    return site::Kind::tail;
  }
  if (name == "donorH") {
    return site::Kind::donorH;
  }
  if (name == "acceptor") {
    return site::Kind::acceptor;
  }
  if (name == "polar") {
    return site::Kind::polar;
  }
  if (name == "apolar") {
    return site::Kind::apolar;
  }
  if (name == "waterO") {
    return site::Kind::waterO;
  }
  if (name == "waterH") {
    return site::Kind::waterH;
  }
  if (name == "solvent") {
    return site::Kind::solvent;
  }
  throw std::invalid_argument("unknown site kind '" + std::string(name) + "'");
}

site::Table parseSiteSpec(std::string_view spec, site::Family family) {
  site::Table table;
  table.family = family;
  std::size_t start = 0;
  while (start <= spec.size()) {
    const std::size_t comma = spec.find(',', start);
    const auto raw =
        spec.substr(start, comma == std::string_view::npos
                               ? std::string_view::npos
                               : comma - start);
    start = (comma == std::string_view::npos) ? spec.size() + 1 : comma + 1;
    const auto token = trimView(raw);
    if (token.empty()) {
      continue;
    }
    const auto eq = token.find('=');
    if (eq == std::string_view::npos) {
      throw std::invalid_argument("site spec token '" + std::string(token) +
                                  "' is not key=value");
    }
    const auto key = trimView(token.substr(0, eq));
    const auto val = trimView(token.substr(eq + 1));
    if (key.empty() || val.empty()) {
      throw std::invalid_argument("site spec token '" + std::string(token) +
                                  "' is not key=value");
    }
    if (key == "family") {
      table.family = site::parseFamily(val);
      continue;
    }
    int typeId = 0;
    try {
      typeId = parseIntegral<int>(key);
    } catch (const ParsingException &) {
      throw std::invalid_argument("site spec type '" + std::string(key) +
                                  "' is not an integer");
    }
    table.typeToKind[typeId] = kindFromName(val);
  }
  return table;
}

void countIonTypes(const Cloud &ions, int &nCation, int &nAnion) {
  nCation = 0;
  nAnion = 0;
  for (const auto &pt : ions.pts) {
    if (pt.type == 1) {
      ++nCation;
    } else if (pt.type == 2) {
      ++nAnion;
    }
  }
}

double meanCageDegree(const Cloud &ions, double cutoff) {
  if (ions.nop == 0) {
    return 0.0;
  }
  const auto nList = nneigh::neighList(cutoff, ions, 1, 2);
  double sum = 0.0;
  int nI = 0;
  const int nRows = static_cast<int>(nList.size());
  for (int i = 0; i < ions.nop; ++i) {
    if (ions.pts[static_cast<std::size_t>(i)].type != 1) {
      continue;
    }
    ++nI;
    if (i >= nRows || nList[static_cast<std::size_t>(i)].empty()) {
      continue;
    }
    sum += static_cast<double>(nList[static_cast<std::size_t>(i)].size() - 1);
  }
  return nI > 0 ? sum / static_cast<double>(nI) : 0.0;
}

int uniqueTypeCount(const Cloud &cloud) {
  std::set<int> types;
  for (const auto &pt : cloud.pts) {
    types.insert(pt.type);
  }
  return static_cast<int>(types.size());
}

int cmdCn(std::ostream &os, Cloud &cloud, double cutoff, int typeI,
          int typeJ) {
  const int typI = typeOf(cloud, typeI);
  const int typJ = typeJ > 0 ? typeJ : typI;
  const auto nList = nneigh::neighListPair(cutoff, cloud, typI, typJ);
  double sum = 0.0;
  int nI = 0;
  const int nRows = static_cast<int>(nList.size());
  for (int i = 0; i < cloud.nop; ++i) {
    if (cloud.pts[static_cast<std::size_t>(i)].type != typI) {
      continue;
    }
    ++nI;
    if (i >= nRows || nList[static_cast<std::size_t>(i)].empty()) {
      continue;
    }
    sum += static_cast<double>(nList[static_cast<std::size_t>(i)].size() - 1);
  }
  const double degree = nI > 0 ? sum / static_cast<double>(nI) : 0.0;
  os << "site-site types " << typI << " " << typJ << " cutoff " << cutoff
     << " degree " << degree << " nI " << nI << "\n";
  return 0;
}

int cmdCnIons(std::ostream &os, Cloud &cloud, double cutoff,
              const site::Table &table) {
  const auto ions = site::ionCloud(cloud, table);
  int nCation = 0;
  int nAnion = 0;
  countIonTypes(ions, nCation, nAnion);
  const double degree = meanCageDegree(ions, cutoff);
  os << "cage ionCloud types 1 2 cutoff " << cutoff << " degree " << degree
     << " nCation " << nCation << " nAnion " << nAnion << "\n";
  return 0;
}

int cmdPairs(std::ostream &os, Cloud &cloud, const site::Table &table) {
  const auto ions = site::ionCloud(cloud, table);
  const auto pairs = nneigh::mutualNearestUnlike(ions, 1, 2);
  int nCation = 0;
  int nAnion = 0;
  countIonTypes(ions, nCation, nAnion);
  os << "contact-pair count " << pairs.size() << " nCation " << nCation
     << " nAnion " << nAnion << "\n";
  return 0;
}

int cmdDomains(std::ostream &os, Cloud &cloud, double cutoff,
               const site::Table &table, site::Kind subset) {
  std::vector<bool> mask(static_cast<std::size_t>(cloud.nop), false);
  for (const int i : site::indicesOf(cloud, table, subset)) {
    if (i >= 0 && i < cloud.nop) {
      mask[static_cast<std::size_t>(i)] = true;
    }
  }
  const auto idxList = nneigh::getNewNeighbourListByIndex(cloud, cutoff);
  std::vector<std::vector<int>> idList(idxList.size());
  for (std::size_t i = 0; i < idxList.size(); ++i) {
    idList[i].reserve(idxList[i].size());
    for (const int idx : idxList[i]) {
      if (idx >= 0 && idx < cloud.nop) {
        idList[i].push_back(cloud.pts[static_cast<std::size_t>(idx)].atomID);
      }
    }
  }
  const auto d = clump::largestDomain(cloud, idList, mask);
  const char *name =
      subset == site::Kind::apolar ? "apolar" : "polar";
  os << "subset " << name << " n " << d.subset << " largest " << d.largest
     << " P " << d.percolation << "\n";
  return 0;
}

int cmdDensityZ(std::ostream &os, Cloud &cloud, int typeI, int bins, int axis) {
  if (bins <= 0) {
    const int a = (axis >= 0 && axis <= 2) ? axis : 2;
    const double span =
        (static_cast<int>(cloud.box.size()) > a)
            ? cloud.box[static_cast<std::size_t>(a)]
            : 0.0;
    bins = std::max(1, static_cast<int>(std::lround(span / 0.1)));
  }
  const auto d = site::densityZ(cloud, typeI, bins, axis);
  os << "# z rho\n";
  for (std::size_t i = 0; i < d.z.size(); ++i) {
    os << d.z[i] << " " << d.rho[i] << "\n";
  }
  return 0;
}

int cmdRdf(std::ostream &os, Cloud &cloud, double rmax, int bins, int typeI,
           int typeJ) {
  if (bins <= 0) {
    bins = std::max(1, static_cast<int>(std::lround(rmax / 0.1)));
  }
  const int typI = typeOf(cloud, typeI);
  const int typJ = typeJ > 0 ? typeJ : typI;
  const auto gr = rdf::partialRdf(cloud, typI, typJ, rmax, bins);
  os << "# r g count\n";
  os << "# types " << typI << " " << typJ << " rmax " << rmax << " bins "
     << bins << " volume " << gr.volume << "\n";
  for (std::size_t i = 0; i < gr.r.size(); ++i) {
    os << gr.r[i] << " " << gr.g[i] << " " << gr.count[i] << "\n";
  }
  return 0;
}

int cmdCn(std::ostream &os, Cloud &cloud, double rmax, int bins, int typeI,
          int typeJ) {
  if (bins <= 0) {
    bins = std::max(1, static_cast<int>(std::lround(rmax / 0.1)));
  }
  const int typI = typeOf(cloud, typeI);
  const int typJ = typeJ > 0 ? typeJ : typI;
  const auto gr = rdf::partialRdf(cloud, typI, typJ, rmax, bins);
  const double rhoJ =
      (gr.volume > 0.0) ? static_cast<double>(gr.nJ) / gr.volume : 0.0;
  const double cn = rdf::coordinationNumber(gr, rmax, rhoJ);
  os << "# site-site\n";
  os << "# types " << typI << " " << typJ << " cutoff " << rmax << " cn "
     << cn << " nI " << gr.nI << " nJ " << gr.nJ << " volume " << gr.volume
     << "\n";
  return 0;
}

int cmdChillPlus(std::ostream &os, Cloud &cloud, double cutoff, int typeI) {
  if (cloud.nop == 0) {
    printCounts(os, cloud);
    return 0;
  }
  const int typ = typeOf(cloud, typeI);
  auto nList = nneigh::neighListO(cutoff, cloud, typ);
  chill::getCorrelPlus(cloud, nList, false);
  chill::getIceTypePlusNoPrint(cloud, nList, false);
  printCounts(os, cloud);
  return 0;
}

int cmdChill(std::ostream &os, Cloud &cloud, double cutoff, int typeI) {
  if (cloud.nop == 0) {
    printCounts(os, cloud);
    return 0;
  }
  const int typ = typeOf(cloud, typeI);
  auto nList = nneigh::neighListO(cutoff, cloud, typ);
  chill::getCorrel(cloud, nList, false);
  chill::getIceTypeNoPrint(cloud, nList, false);
  printCounts(os, cloud);
  return 0;
}

int cmdCages(std::ostream &os, Cloud &cloud, double cutoff, int typeI, int k,
             const std::string &graphName) {
  if (cloud.nop == 0) {
    os << colorizer.heading("nop") << " 0 "
       << colorizer.longOption("graph") << " " << graphName << " "
       << iceColor("hexagonal") << " 0 " << iceColor("cubic") << " 0 "
       << iceColor("water") << " 0\n";
    return 0;
  }
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
    auto graphs = nneigh::kNearestNeighbourPair(cloud, k, cand, typ);
    const auto &mutual = graphs.first;
    const auto &uni = graphs.second;
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
  os << colorizer.heading("nop") << " " << cloud.nop << " "
     << colorizer.longOption("graph") << " " << graphName << " "
     << iceColor("hexagonal") << " " << ih << " " << iceColor("cubic") << " "
     << ic << " " << iceColor("water") << " " << water << "\n";
  return 0;
}

int cmdHbonds(std::ostream &os, Cloud &yCloud, Cloud &hCloud, double cutoff,
              int typeI, double distCutoff, double angleCutoff,
              bool allDonors) {
  if (yCloud.nop == 0) {
    os << colorizer.heading("nop") << " 0 "
       << colorizer.longOption("hbonds") << " 0\n";
    return 0;
  }
  const int typ = typeOf(yCloud, typeI);
  auto nList = nneigh::neighListO(cutoff, yCloud, typ);
  std::vector<std::vector<int>> net;
  if (allDonors) {
    std::vector<int> donorHs;
    donorHs.reserve(static_cast<std::size_t>(hCloud.nop));
    for (int i = 0; i < hCloud.nop; ++i) {
      donorHs.push_back(i);
    }
    net = bond::populateHbondsFromDonors(yCloud, hCloud, nList, donorHs,
                                         distCutoff, angleCutoff);
  } else {
    net = bond::populateHbondsWithInputClouds(yCloud, hCloud, nList, distCutoff,
                                              angleCutoff);
  }
  int edges = 0;
  for (const auto &row : net) {
    if (row.size() > 1) {
      edges += static_cast<int>(row.size()) - 1;
    }
  }
  os << colorizer.heading("nop") << " " << yCloud.nop << " "
     << colorizer.longOption("hbonds") << " " << (edges / 2) << "\n";
  return 0;
}

} // namespace

int main(int argc, char *argv[]) {
  colorizer = colorizerForFile(environmentColorStatus(), stdout);

  auto cfg = seams::cfg::load(seams::cfg::pathFromArgv(argc, argv));
  int frame = cfg.frame;
  int last = cfg.last;
  int jobs = cfg.jobs;
  int typeI = cfg.type;
  double cutoff = cfg.cutoff;
  int k = cfg.k;
  std::string graph = cfg.graph;
  site::Family family = cfg.family;
  bool familyFlag = false;
  std::string familyNameFlag;
  bool printConfig = false;
  bool strictInput = false;
  int htype = 1;
  double hdist = 2.42;
  double hangle = 30.0;
  bool allDonors = false;
  std::string cmd;
  std::string file;
  std::string typesFlag;
  std::string siteSpec;
  bool ionsFlag = false;
  std::string subsetFlag;
  int rdfTypeI = 0;
  int rdfTypeJ = 0;
  bool typesSet = false;
  int bins = 0;
  int densAxis = 2;
  std::string axisFlag;
  std::string outputFormat = "text";

  const char *progname = (argc ? argv[0] : "seams");
  Parser parser;

  parser.add(Option("--help", "-h")
                 .help("show this help message and exit")
                 .handler([&]() {
                   std::cout << colorizer.heading("d-SEAMS") << "\n\n";
                   std::cout << parser.formatHelp(progname, colorizer);
                   std::exit(EXIT_SUCCESS);
                 }));

  parser.add(Option("--version", "-v")
                 .help("Print version information")
                 .handler([&]() {
#ifdef SEAMS_VERSION
                   std::cout << colorizer.progName("seams") << " "
                             << SEAMS_VERSION << "\n";
#else
                   std::cout << colorizer.progName("seams") << "\n";
#endif
                   std::exit(EXIT_SUCCESS);
                 }));

  parser.add(Option("--features")
                 .help("Print compile-time backends")
                 .handler([&]() {
                   printFeatures(std::cout);
                   std::exit(EXIT_SUCCESS);
                 }));

  parser.add(Option("--format")
                 .argName("KIND")
                 .help("Output format: text (default) or json")
                 .handler([&](std::string_view value) {
                   outputFormat = std::string(value);
                 }));

  parser.add(Option("--strict-input")
                 .help("Fail when a requested frame produces no atoms")
                 .handler([&]() { strictInput = true; }));

  parser.add(Option("--frame", "-f")
                 .argName("N")
                 .help("First frame (1-based)")
                 .handler([&](std::string_view value) {
                   frame = parseIntegral<int>(value);
                 }));

  parser.add(Option("--last")
                 .argName("N")
                 .help("Last frame (inclusive). Omit for a single --frame")
                 .handler([&](std::string_view value) {
                   last = parseIntegral<int>(value);
                 }));

  parser.add(Option("--jobs", "-j")
                 .argName("N")
                 .help("Parallel frame workers (OpenMP). 1 is serial")
                 .handler([&](std::string_view value) {
                   jobs = parseIntegral<int>(value);
                 }));

  parser.add(Option("--type", "-t")
                 .argName("I")
                 .help("Atom type (0 guesses oxygen then type 1; "
                       "density-z: 0 is every atom)")
                 .handler([&](std::string_view value) {
                   typeI = parseIntegral<int>(value);
                 }));

  parser.add(Option("--cutoff", "-c")
                 .argName("ANGSTROM")
                 .help("Neighbour cutoff (rdf/cn rmax; hbonds heavy-atom nList)")
                 .handler([&](std::string_view value) {
                   cutoff = parseFloatingPoint<double>(value);
                 }));

  parser.add(Option("--types")
                 .argName("I,J")
                 .help("Pair types for site-site rdf/cn (default I=J=--type)")
                 .handler([&](std::string_view value) {
                   typesFlag = std::string(value);
                 }));

  parser.add(Option("--family")
                 .argName("NAME")
                 .help("Material family (default waterIce): waterIce, "
                       "ionicLiquid, moltenSalt, des, electrolyte, "
                       "confinedIL, confinedWater, networkFormer")
                 .handler([&](std::string_view value) {
                   familyNameFlag = std::string(value);
                   familyFlag = true;
                 }));

  parser.add(Option("--site")
                 .argName("SPEC")
                 .help("Type-to-kind map for ionCloud (1=cationHead,2=anion)")
                 .handler([&](std::string_view value) {
                   siteSpec = std::string(value);
                 }));

  parser.add(Option("--ions")
                 .help("cn on ionCloud cage degree (needs --site)")
                 .handler([&]() { ionsFlag = true; }));

  parser.add(Option("--subset")
                 .argName("KIND")
                 .help("domains subset: polar | apolar")
                 .handler([&](std::string_view value) {
                   subsetFlag = std::string(value);
                 }));

  parser.add(Option("--bins")
                 .argName("N")
                 .help("Histogram bins (rdf: rmax/0.1; density-z: L/0.1)")
                 .handler([&](std::string_view value) {
                   bins = parseIntegral<int>(value);
                 }));

  parser.add(Option("--axis")
                 .argName("XYZ")
                 .help("Histogram axis for density-z (x|y|z, default z)")
                 .handler([&](std::string_view value) {
                   axisFlag = std::string(value);
                 }));

  parser.add(Option("-k")
                 .argName("N")
                 .help("k for knn / seeded cages")
                 .handler([&](std::string_view value) {
                   k = parseIntegral<int>(value);
                 }));

  parser.add(Option("--graph")
                 .argName("KIND")
                 .help("Bond graph for cages: cutoff | knn | knn-union | seeded")
                 .handler([&](std::string_view value) { graph = value; }));

  parser.add(Option("--htype")
                 .argName("I")
                 .help("Hydrogen atom type for hbonds")
                 .handler([&](std::string_view value) {
                   htype = parseIntegral<int>(value);
                 }));

  parser.add(Option("--hdist")
                 .argName("ANGSTROM")
                 .help("Acceptor-H distance cutoff for hbonds")
                 .handler([&](std::string_view value) {
                   hdist = parseFloatingPoint<double>(value);
                 }));

  parser.add(Option("--hangle")
                 .argName("DEG")
                 .help("O-O-H angle cutoff for hbonds (acceptor-centered)")
                 .handler([&](std::string_view value) {
                   hangle = parseFloatingPoint<double>(value);
                 }));

  parser.add(Option("--donors")
                 .help("Use every hydrogen as a donor candidate")
                 .handler([&]() { allDonors = true; }));

  parser.add(Option("--config")
                 .argName("FILE")
                 .help("Dotenv file of SEAMS_* / LINKCELL_* knobs "
                       "(or set SEAMS_CONFIG). Env wins over the file")
                 .handler([&](std::string_view) {}));

  parser.add(Option("--print-config")
                 .help("Print the resolved runtime knobs and exit")
                 .handler([&]() { printConfig = true; }));

  parser.add(Option("--tpp")
                 .argName("N")
                 .help("Threads per particle for the device k-NN "
                       "(LINKCELL_TPP). 0 is the occupancy picker")
                 .handler([&](std::string_view value) {
                   cfg.tpp = parseIntegral<int>(value);
                 }));

  parser.add(Option("--block")
                 .argName("N")
                 .help("CUDA block size (LINKCELL_BLOCK / SEAMS_BLOCK)")
                 .handler([&](std::string_view value) {
                   cfg.block = parseIntegral<int>(value);
                 }));

  parser.add(Option("--resident")
                 .argName("FRAC")
                 .help("Fraction of free device memory for a resident batch "
                       "(SEAMS_RESIDENT, default 0.80)")
                 .handler([&](std::string_view value) {
                   cfg.resident = parseFloatingPoint<double>(value);
                 }));

  parser.add(Option("--cell")
                 .argName("ANGSTROM")
                 .help("Link-cell hint so NPT frames share a grid "
                       "(SEAMS_CELL, default 3)")
                 .handler([&](std::string_view value) {
                   cfg.cell = parseFloatingPoint<double>(value);
                 }));

  parser.add(Positional("command")
                 .help("read | chill | chill-plus | cages | rdf | cn | hbonds | "
                       "pairs | density-z | domains")
                 .occurs(zeroOrOneTime)
                 .handler([&](std::string_view value) { cmd = value; }));

  parser.add(Positional("file")
                 .help("Trajectory file")
                 .occurs(zeroOrOneTime)
                 .handler([&](std::string_view value) { file = value; }));

  try {
    parser.parse(argc, argv);
  } catch (const ParsingException &ex) {
    std::cerr << colorizer.error(ex.message()) << "\n";
    std::cerr << colorizer.warning(parser.formatUsage(progname, colorizer))
              << "\n";
    return 2;
  }

  if (outputFormat != "text" && outputFormat != "json") {
    std::cerr << colorizer.error("bad --format (want text|json)") << "\n";
    return 2;
  }

  cfg.frame = frame;
  cfg.last = last;
  cfg.jobs = jobs;
  cfg.type = typeI;
  cfg.cutoff = cutoff;
  cfg.k = k;
  cfg.graph = graph;
  if (familyFlag) {
    try {
      family = site::parseFamily(familyNameFlag);
    } catch (const std::invalid_argument &e) {
      std::cerr << colorizer.error(e.what()) << "\n";
      std::cerr << colorizer.warning(parser.formatUsage(progname, colorizer))
                << "\n";
      return 2;
    }
  }
  cfg.family = family;
  seams::cfg::exportEnviron(cfg);

  if (printConfig) {
    seams::cfg::dump(cfg, std::cout);
    return 0;
  }

  if (!typesFlag.empty()) {
    if (!parseTypePair(typesFlag, rdfTypeI, rdfTypeJ)) {
      std::cerr << colorizer.error("bad --types (want I,J)") << "\n";
      std::cerr << colorizer.warning(parser.formatUsage(progname, colorizer))
                << "\n";
      return 2;
    }
    typesSet = true;
  } else {
    rdfTypeI = typeI;
    rdfTypeJ = typeI;
  }
  if (bins < 0) {
    std::cerr << colorizer.error("bad --bins (want N > 0)") << "\n";
    return 2;
  }
  if (!axisFlag.empty()) {
    const auto axisView = trimView(axisFlag);
    if (axisView == "x" || axisView == "X" || axisView == "0") {
      densAxis = 0;
    } else if (axisView == "y" || axisView == "Y" || axisView == "1") {
      densAxis = 1;
    } else if (axisView == "z" || axisView == "Z" || axisView == "2") {
      densAxis = 2;
    } else {
      std::cerr << colorizer.error("bad --axis (want x|y|z)") << "\n";
      std::cerr << colorizer.warning(parser.formatUsage(progname, colorizer))
                << "\n";
      return 2;
    }
  }

  if (cmd.empty() || file.empty()) {
    std::cerr << colorizer.error(
                     "A command and a trajectory file are required")
              << "\n";
    std::cerr << colorizer.warning(parser.formatUsage(progname, colorizer))
              << "\n";
    return 2;
  }

  const bool iceCmd = cmd == "chill" || cmd == "chill-plus" ||
                      cmd == "chill_plus" || cmd == "cages";
  if (iceCmd && !site::iceScoreAllowed(family)) {
    std::cerr << colorizer.error(site::refuseIceScore(family)) << "\n";
    return 2;
  }

  site::Table siteTable;
  siteTable.family = family;
  if (!siteSpec.empty()) {
    try {
      siteTable = parseSiteSpec(siteSpec, family);
    } catch (const std::exception &e) {
      std::cerr << colorizer.error(e.what()) << "\n";
      return 2;
    }
  }
  if (cmd == "pairs" && siteSpec.empty()) {
    std::cerr << colorizer.error("pairs needs --site") << "\n";
    return 2;
  }
  if (cmd == "cn" && ionsFlag && siteSpec.empty()) {
    std::cerr << colorizer.error("cn --ions needs --site") << "\n";
    return 2;
  }
  if (cmd == "domains" && siteSpec.empty()) {
    std::cerr << colorizer.error("domains needs --site") << "\n";
    return 2;
  }
  site::Kind domainKind = site::Kind::polar;
  if (cmd == "domains") {
    if (subsetFlag == "apolar" || subsetFlag == "tail") {
      domainKind = site::Kind::apolar;
    } else if (subsetFlag.empty() || subsetFlag == "polar") {
      domainKind = site::Kind::polar;
    } else {
      std::cerr << colorizer.error("bad --subset (want polar|apolar)") << "\n";
      return 2;
    }
  }

  auto runOne = [&](std::ostream &os, Cloud &cloud) {
    if (strictInput && cloud.nop == 0) {
      os << "input frame " << cloud.currentFrame << " contains no atoms\n";
      return 2;
    }
    if (!familyFlag && family == site::Family::waterIce &&
        uniqueTypeCount(cloud) > 2) {
      std::cerr << colorizer.warning(
                       "hint: more than two LAMMPS types; set --family "
                       "(default waterIce)")
                << "\n";
    }
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
        os << colorizer.error(e.what()) << "\n";
        return 2;
      }
    }
    if (cmd == "rdf") {
      return cmdRdf(os, cloud, cutoff, bins, rdfTypeI, rdfTypeJ);
    }
    if (cmd == "cn") {
      if (ionsFlag) {
        return cmdCnIons(os, cloud, cutoff, siteTable);
      }
      return cmdCn(os, cloud, cutoff, bins, rdfTypeI, rdfTypeJ);
    }
    if (cmd == "pairs") {
      return cmdPairs(os, cloud, siteTable);
    }
    if (cmd == "density-z") {
      return cmdDensityZ(os, cloud, typeI, bins, densAxis);
    }
    if (cmd == "domains") {
      return cmdDomains(os, cloud, cutoff, siteTable, domainKind);
    }
    if (cmd == "hbonds") {
      Cloud hCloud = load(file, cloud.currentFrame, htype);
      return cmdHbonds(os, cloud, hCloud, cutoff, typeI, hdist, hangle,
                       allDonors);
    }
    os << colorizer.error("unknown command: ") << cmd << "\n";
    return 2;
  };

  const bool loadAll =
      cmd == "pairs" || cmd == "hbonds" || (cmd == "cn" && ionsFlag) ||
      ((cmd == "rdf" || cmd == "cn") && typesSet && rdfTypeI != rdfTypeJ) ||
      (cmd == "density-z" && typeI <= 0) || cmd == "domains";
  const int loadType = loadAll ? -1 : typeI;

  if (last <= 0 || last == frame) {
    Cloud cloud = load(file, frame, loadType);
    std::ostringstream line;
    const int rc = runOne(line, cloud);
    emitResult(std::cout, outputFormat, cmd, cloud.currentFrame, rc,
               line.str());
    return rc;
  }

  std::mutex outMu;
  const int outputCount = std::max(0, last - frame + 1);
  std::vector<std::string> outputLines(static_cast<std::size_t>(outputCount));
  std::vector<int> outputStatus(static_cast<std::size_t>(outputCount), 0);
  std::vector<bool> outputSeen(static_cast<std::size_t>(outputCount), false);
  const int typeFilter = loadAll ? 0 : (typeI > 0 ? typeI : 0);
  auto processFrame = [&](int fr, Cloud &cloud) {
        if (typeFilter <= 0 && cloud.nop == 0) {
          cloud = load(file, cloud.currentFrame, loadType);
        }
        std::ostringstream line;
        const int one = runOne(line, cloud);
        {
          std::lock_guard<std::mutex> lock(outMu);
          const int index = fr - frame;
          if (index >= 0 && index < outputCount) {
            outputLines[static_cast<std::size_t>(index)] = line.str();
            outputStatus[static_cast<std::size_t>(index)] = one;
            outputSeen[static_cast<std::size_t>(index)] = true;
          }
        }
      };
  if (!forEachInputFrame(file, frame, last, typeFilter, jobs, processFrame)) {
    std::cerr << colorizer.error("frame ranges are unsupported for this input format")
              << "\n";
    return 2;
  }
  int rc = 0;
  for (std::size_t i = 0; i < outputLines.size(); ++i) {
    if (!outputSeen[i]) {
      continue;
    }
    const int frameNumber = frame + static_cast<int>(i);
    emitResult(std::cout, outputFormat, cmd, frameNumber, outputStatus[i],
               outputLines[i]);
    if (rc == 0 && outputStatus[i] != 0) {
      rc = outputStatus[i];
    }
  }
  return rc;
}
