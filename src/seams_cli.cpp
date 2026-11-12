//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <argum.h>
#include <bond.hpp>
#include <bop.hpp>
#include <chill_plus_scheme.hpp>
#include <cage_affiliation.hpp>
#include <cage_enum.hpp>
#include <franzblau.hpp>
#include <local_symmetry.hpp>
#include <mol_sys.hpp>
#include <generic.hpp>
#include <neighbours.hpp>
#include <order_parameter.hpp>
#include <cluster.hpp>
#include <density.hpp>
#include <edge_set.hpp>
#include <rdf.hpp>
#include <seams_config.hpp>
#include <seams_input.hpp>
#include <site.hpp>
#include <topo_bulk.hpp>
#include <water_select.hpp>
#include <topo_fingerprint.hpp>
#include <tum_offload.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#ifdef SEAMS_HAS_MPI
#include <mpi.h>
#include <unistd.h>

#include <cstring>
#include <limits>
#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif
#endif

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
  if (ext == "con") {
#ifdef SEAMS_HAS_READCON
    // readCon prints an open error and returns an empty cloud. A missing
    // path has to fail the command, or the CLI reports success.
    if (!gen::file_exists(path)) {
      throw std::runtime_error("cannot read " + path);
    }
    return sinp::readCon(path, frame, cloud);
#else
    throw std::runtime_error("readcon-core is not in this build; CON input needs it");
#endif
  }
  if (ext == "pdb" || ext == "gro" || ext == "dcd") {
#ifdef SEAMS_HAS_CHEMFILES
    const int filter = typeI < 0 ? -1 : (typeI > 0 ? typeI : -1);
    return sinp::readChemfiles(path, frame, cloud, filter);
#else
    throw std::runtime_error("chemfiles is not in this build; " + ext + " input needs it");
#endif
  }
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

//! How a frame range is read: an indexed LAMMPS dump, one XYZ frame, or a
//! reader that loads frame by frame until one is missing
enum class FrameSource { dump, xyz, serial };

FrameSource frameSource(const std::string &path) {
  const auto dot = path.find_last_of('.');
  const std::string ext = dot == std::string::npos ? "" : path.substr(dot + 1);
  if (ext == "xyz") {
    return FrameSource::xyz;
  }
  if (ext == "con" || ext == "pdb" || ext == "gro" || ext == "dcd") {
    return FrameSource::serial;
  }
  return FrameSource::dump;
}

template <typename Fn>
bool forEachInputFrame(const std::string &path, int first, int last,
                       int typeFilter, int nThreads, Fn &&fn, int part = 0,
                       int parts = 1) {
  const FrameSource source = frameSource(path);
  if (source == FrameSource::xyz) {
    // last is an int, so last <= 0 || last >= 1 holds for every value.
    if (first <= 1) {
      Cloud cloud = load(path, 1, typeFilter);
      fn(1, cloud);
    }
    return true;
  }

  if (source == FrameSource::serial) {
    const int stop = last;
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
                           std::forward<Fn>(fn), nThreads, part, parts);
  return true;
}

#ifdef SEAMS_HAS_MPI
//! Ranks of an MPI launch, sharing a dump's frames in serpentine rounds.
//! A bare run, or SEAMS_MPI=0, stays out of MPI: a singleton start costs
//! about a second, and a seams started inside an MPI rank inherits its
//! launcher variables.
struct FrameRanks {
  bool active = false;
  bool sameCommand = true;
  int rank = 0;
  int size = 1;

  FrameRanks(int *argc, char ***argv) {
    const char *opt = std::getenv("SEAMS_MPI");
    if (opt != nullptr && std::string_view(opt) == "0") {
      return;
    }
    int launched = 0;
    for (const char *name :
         {"OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "MV2_COMM_WORLD_SIZE"}) {
      if (const char *value = std::getenv(name)) {
        launched = std::max(launched, std::atoi(value));
      }
    }
    active = launched > 0 || std::getenv("PMIX_RANK") != nullptr ||
             std::getenv("PMI_RANK") != nullptr;
    if (!active) {
      return;
    }
    // Ranks sharing a node would each start a thread per core; rayon reads
    // its count once, and MPI's own threads must not race a setenv
    if (launched > 1) {
      setenv("RAYON_NUM_THREADS", "1", 0);
    }
    int provided = 0;
    MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#ifdef SEAMS_HAS_OPENMP
    if (size > 1 && std::getenv("OMP_NUM_THREADS") == nullptr) {
      omp_set_num_threads(1);
    }
#endif
    // Every rank must analyse one command on one input, or rank 0 prints
    // another rank's frames as its own: hash the arguments, the directory,
    // and every SEAMS_ variable, which by now hold the config file too
    std::vector<std::string> said(*argv + 1, *argv + *argc);
    char cwd[4096] = {};
    said.emplace_back(getcwd(cwd, sizeof cwd) != nullptr ? cwd : "");
    std::vector<std::string> settings;
    for (char **e = environ; *e != nullptr; e++) {
      if (std::string_view(*e).rfind("SEAMS_", 0) == 0) {
        settings.emplace_back(*e);
      }
    }
    std::sort(settings.begin(), settings.end());
    said.insert(said.end(), settings.begin(), settings.end());
    unsigned long long hash = 14695981039346656037ULL;
    for (const std::string &s : said) {
      for (const char c : s) {
        hash = (hash ^ static_cast<unsigned char>(c)) * 1099511628211ULL;
      }
      hash = (hash ^ 0xffULL) * 1099511628211ULL;
    }
    unsigned long long lowest = 0;
    unsigned long long highest = 0;
    MPI_Allreduce(&hash, &lowest, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&hash, &highest, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX,
                  MPI_COMM_WORLD);
    sameCommand = lowest == highest;
    // A rank analyses whole frames, so steinhardtQl keeps its atoms on it
    chill::setSteinhardtAtomSplit(false);
    if (rank != 0) {
      std::cout.setstate(std::ios::failbit);
    }
  }
  ~FrameRanks() {
    if (active) {
      MPI_Finalize();
    }
  }
  FrameRanks(const FrameRanks &) = delete;
  FrameRanks &operator=(const FrameRanks &) = delete;
};

//! Every rank's (frame index, status, text) records, unpacked on rank 0.
//! False on every rank when they would pass MPI's int counts, and on rank 0
//! when a record does not unpack.
bool gatherFrameResults(const FrameRanks &ranks,
                        std::vector<std::string> &lines,
                        std::vector<int> &status, std::vector<bool> &seen) {
  std::string pack;
  for (std::size_t i = 0; i < lines.size(); i++) {
    if (seen[i]) {
      const int head[3] = {static_cast<int>(i), status[i],
                           static_cast<int>(lines[i].size())};
      pack.append(reinterpret_cast<const char *>(head), sizeof head);
      pack += lines[i];
    }
  }
  const long long bytes = static_cast<long long>(pack.size());
  long long total = 0;
  MPI_Allreduce(&bytes, &total, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  if (total > std::numeric_limits<int>::max()) {
    return false;
  }
  const int mine = static_cast<int>(bytes);
  std::vector<int> counts(static_cast<std::size_t>(ranks.size), 0);
  MPI_Gather(&mine, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
  std::vector<int> displs(counts.size(), 0);
  std::string all;
  if (ranks.rank == 0) {
    for (std::size_t r = 1; r < counts.size(); r++) {
      displs[r] = displs[r - 1] + counts[r - 1];
    }
    all.resize(static_cast<std::size_t>(total));
  }
  MPI_Gatherv(pack.data(), mine, MPI_BYTE, all.data(), counts.data(),
              displs.data(), MPI_BYTE, 0, MPI_COMM_WORLD);
  for (std::size_t at = 0; at < all.size();) {
    int head[3];
    if (all.size() - at < sizeof head) {
      return false;
    }
    std::memcpy(head, all.data() + at, sizeof head);
    at += sizeof head;
    const auto i = static_cast<std::size_t>(head[0]);
    const auto len = static_cast<std::size_t>(head[2]);
    if (head[0] < 0 || head[2] < 0 || i >= lines.size() ||
        len > all.size() - at) {
      return false;
    }
    lines[i] = all.substr(at, len);
    status[i] = head[1];
    seen[i] = true;
    at += len;
  }
  return true;
}
#endif

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

std::string exceptionText() {
  try {
    throw;
  } catch (const ParsingException &ex) {
    return std::string(ex.message());
  } catch (const std::exception &ex) {
    return ex.what();
  } catch (...) {
    return "unknown error";
  }
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

bool isWaterType(int type, const std::vector<int> &waterTypes) {
  if (waterTypes.empty()) {
    return true;
  }
  return std::find(waterTypes.begin(), waterTypes.end(), type) !=
         waterTypes.end();
}

int waterNop(const Cloud &cloud, const std::vector<int> &waterTypes) {
  if (waterTypes.empty()) {
    return cloud.nop;
  }
  int n = 0;
  for (const auto &pt : cloud.pts) {
    if (isWaterType(pt.type, waterTypes)) {
      ++n;
    }
  }
  return n;
}

std::vector<std::vector<int>>
waterGraph(Cloud &cloud, double cutoff, int typ,
           const std::vector<int> &waterTypes, int k = 4, bool mutual = true) {
  const double cand = cutoff + 1.5;
  if (!waterTypes.empty()) {
    return nneigh::kNearestNeighbourList(cloud, k, cand, waterTypes, mutual);
  }
  return nneigh::kNearestNeighbourList(cloud, k, cand, typ, mutual);
}

std::string iceColor(std::string_view name) {
  if (name == "cubic") {
    return colorizer.longOption(name);
  }
  if (name == "hexagonal") {
    return colorizer.shortOption(name);
  }
  if (name == "interfacial") {
    return colorizer.warning(name);
  }
  if (name == "clathrate" || name == "interClathrate") {
    return colorizer.progName(name);
  }
  if (name == "water" || name == "both") {
    return std::string(name);
  }
  return colorizer.error(name);
}

void printCounts(std::ostream &os, const Cloud &cloud,
                 const std::vector<int> &waterTypes = {}) {
  std::map<std::string, int> hist;
  int scored = 0;
  for (const auto &pt : cloud.pts) {
    if (!isWaterType(pt.type, waterTypes)) {
      continue;
    }
    hist[iceName(pt.iceType)]++;
    ++scored;
  }
  os << colorizer.heading("nop") << " " << scored;
  for (const auto &[name, n] : hist) {
    os << " " << iceColor(name) << " " << n;
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
#ifdef SEAMS_HAS_READCON_DB
  line("readcon-db", true);
#else
  line("readcon-db", false);
#endif
#ifdef SEAMS_HAS_IRA
  line("IRA/SOFI", true);
#else
  line("IRA/SOFI", false);
#endif
#ifdef SEAMS_HAS_SPGLIB
  line("spglib", true);
#else
  line("spglib", false);
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

void writeOverlay(std::ostream &os, const char *label, const symm::Overlay &hit) {
  os << label << " " << hit.name << " rmsd " << hit.rmsd << " hausdorff "
     << hit.hausdorff << " backend " << hit.backend << "\n";
  os << "rotation";
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      os << " " << hit.rotation(i, j);
    }
  }
  os << "\ntranslation " << hit.translation[0] << " " << hit.translation[1]
     << " " << hit.translation[2] << "\nperm";
  for (int p : hit.permutation) {
    os << " " << p;
  }
  os << "\n";
}

int cmdSymmetry(std::ostream &os, Cloud &cloud) {
  if (cloud.nop < 1 || cloud.box.size() < 3) {
    os << "symmetry needs a frame and a box\n";
    return 2;
  }
  Eigen::MatrixXd xyz(cloud.nop, 3);
  for (int i = 0; i < cloud.nop; i++) {
    const auto &p = cloud.pts[static_cast<size_t>(i)];
    xyz(i, 0) = p.x;
    xyz(i, 1) = p.y;
    xyz(i, 2) = p.z;
  }
  Eigen::Matrix3d vectors = Eigen::Matrix3d::Zero();
  vectors(0, 0) = cloud.box[0];
  vectors(1, 1) = cloud.box[1];
  vectors(2, 2) = cloud.box[2];
  os << "symmetry backend " << symm::matcher() << "\n";
  const symm::Space space = symm::spaceGroup(vectors, xyz, 1e-3);
  if (space.status == 0) {
    os << "space-group " << space.number << " " << space.symbol << "\n";
  } else if (space.status == 2) {
    os << "space-group off\n";
  } else {
    os << "space-group none\n";
  }
  const symm::Overlay whole = symm::best(xyz);
  if (whole.status == 0) {
    writeOverlay(os, "whole", whole);
    const symm::Group group = symm::pointGroup(xyz, 0.05);
    os << "point-group " << (group.status == 0 ? group.symbol : "off") << " "
       << group.nOperations << "\n";
  }
  std::map<std::string, int> tally;
  int unmatched = 0;
  bool wroteShell = false;
  std::vector<int> kind;
  symm::classifyIce(xyz, vectors, kind);
  for (int i = 0; i < cloud.nop; i++) {
    if (kind[static_cast<size_t>(i)] < 0) {
      unmatched++;
      continue;
    }
    const char *label =
        symm::name(static_cast<symm::Kind>(kind[static_cast<size_t>(i)]));
    tally[label] += 1;
    if (!wroteShell) {
      const Eigen::MatrixXd shell = symm::iceShell(xyz, vectors, i);
      const symm::Overlay hit = symm::best(shell);
      if (hit.status == 0) {
        writeOverlay(os, "shell", hit);
        wroteShell = true;
      }
    }
  }
  os << "shells";
  for (const auto &item : tally) {
    os << " " << item.first << " " << item.second;
  }
  os << " unmatched " << unmatched << "\n";
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
    } else {
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
    mask[static_cast<std::size_t>(i)] = true;
  }
  const auto idxList = nneigh::getNewNeighbourListByIndex(cloud, cutoff);
  std::vector<std::vector<int>> idList(idxList.size());
  for (std::size_t i = 0; i < idxList.size(); ++i) {
    idList[i].reserve(idxList[i].size());
    for (const int idx : idxList[i]) {
      idList[i].push_back(cloud.pts[static_cast<std::size_t>(idx)].atomID);
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

//! The graphs one command builds on a frame, each with its own carried state
enum class GraphSlot { single, strict, permissive };

/** Rings and cage affiliation carried from one frame to the next.
 *  A frame walk on one worker visits frames in order through a single
 *  instance, so RingUpdater and AffiliationUpdater redo only the part of
 *  each frame near a changed bond. Their results equal ringNetwork and
 *  cageAffiliation on the same frame. Single frames, parallel workers and
 *  --rebuild use an instance with carry off, which recomputes every call.
 */
class FrameCarry {
public:
  using Rows = std::vector<std::vector<int>>;

  explicit FrameCarry(bool carry) : carry_(carry) {}

  const Rows &rings(GraphSlot slot, const Rows &idx, int depth) {
    const int key = static_cast<int>(slot);
    if (!carry_) {
      return scratchRings_[key] = primitive::ringNetwork(idx, depth);
    }
    auto it = ringUpdaters_.find({key, depth});
    if (it == ringUpdaters_.end()) {
      it = ringUpdaters_.emplace(std::make_pair(key, depth),
                                 primitive::RingUpdater(depth))
               .first;
    }
    return it->second.update(idx);
  }

  const ring::CageAffiliation &affiliation(GraphSlot slot, const Rows &six,
                                           const Rows &idx) {
    const int key = static_cast<int>(slot);
    if (!carry_) {
      return scratchAffiliation_[key] = ring::cageAffiliation(six, idx);
    }
    return affiliationUpdaters_[key].update(six, idx);
  }

private:
  bool carry_;
  std::map<std::pair<int, int>, primitive::RingUpdater> ringUpdaters_;
  std::map<int, ring::AffiliationUpdater> affiliationUpdaters_;
  std::map<int, Rows> scratchRings_;
  std::map<int, ring::CageAffiliation> scratchAffiliation_;
};

std::vector<std::vector<int>> sixRings(const std::vector<std::vector<int>> &rings) {
  std::vector<std::vector<int>> six;
  for (const auto &r : rings) {
    if (r.size() == 6) {
      six.push_back(r);
    }
  }
  return six;
}

// Label-independent topology keys of the bonded graph: one key per atom
// for the rooted neighbourhood within `hops` bonds, the histogram of keys,
// the primitive ring census, and a frame key over all of it.
// One LAMMPS dump frame with an extra per-atom column, so OVITO, VMD or a
// script can colour a trajectory by what the engine decided. `values` is
// indexed like `cloud.pts`; missing entries write 0.
void writePerAtomDump(const std::string &path, const Cloud &cloud, const std::string &column,
                      const std::vector<int> &values) {
  std::ofstream out(path, std::ios::app);
  if (!out) {
    throw std::runtime_error("cannot write per-atom dump " + path);
  }
  out << "ITEM: TIMESTEP\n" << cloud.currentFrame << "\n";
  out << "ITEM: NUMBER OF ATOMS\n" << cloud.nop << "\n";
  gen::writeDumpBoxBounds(out, cloud);
  out << "ITEM: ATOMS id type x y z " << column << "\n";
  for (int i = 0; i < cloud.nop; ++i) {
    const auto &pt = cloud.pts[static_cast<std::size_t>(i)];
    const int v = values[static_cast<std::size_t>(i)];
    out << pt.atomID << ' ' << pt.type << ' ' << pt.x << ' ' << pt.y << ' ' << pt.z << ' ' << v
        << '\n';
  }
}

int cmdFingerprint(std::ostream &os, FrameCarry &carry, Cloud &cloud,
                   double cutoff, int typeI, int k,
                   const std::string &graphName, int hops, bool colourTypes,
                   const std::string &libraryPath, const std::string &emitLabel,
                   const std::string &perAtomPath = {}, bool byAtom = false) {
  if (cloud.nop == 0) {
    os << colorizer.heading("nop") << " 0\n";
    return 0;
  }
  const int typ = typeOf(cloud, typeI);
  const double cand = cutoff + 1.5;
  // the seeded assignment is two graphs; a fingerprint wants one, so the
  // cages default falls back to the cutoff graph here
  const std::string name =
      (graphName.empty() || graphName == "seeded") ? std::string("cutoff") : graphName;
  const auto graph = nneigh::bondGraphFromName(name);
  std::vector<std::vector<int>> nList;
  if (graph == nneigh::BondGraph::Cutoff) {
    nList = nneigh::neighListO(cutoff, cloud, typ);
  } else {
    const bool mutual = graph == nneigh::BondGraph::KnnMutual;
    nList = nneigh::kNearestNeighbourList(cloud, k, cand, typ, mutual);
  }
  const auto rows = nneigh::neighbourListByIndex(cloud, nList);
  std::vector<int> colours;
  if (colourTypes) {
    colours.reserve(static_cast<std::size_t>(cloud.nop));
    for (int i = 0; i < cloud.nop; ++i) {
      colours.push_back(cloud.pts[static_cast<std::size_t>(i)].type);
    }
    colours.resize(rows.size(), 0);
  }
  const auto fp = topo::fingerprint(rows, hops, 7, colours);
  if (!emitLabel.empty()) {
    topo::KeyLibrary lib;
    topo::addToLibrary(lib, fp, emitLabel);
    os << topo::writeLibrary(lib);
    return 0;
  }
  if (!libraryPath.empty()) {
    // several libraries, comma separated, at different depths: the
    // deepest that knows an atom names it
    std::vector<topo::KeyLibrary> libs;
    std::size_t start = 0;
    while (start <= libraryPath.size()) {
      std::size_t comma = libraryPath.find(',', start);
      if (comma == std::string::npos) {
        comma = libraryPath.size();
      }
      const std::string path = libraryPath.substr(start, comma - start);
      start = comma + 1;
      if (path.empty()) {
        continue;
      }
      std::ifstream in(path);
      if (!in) {
        throw std::runtime_error("cannot read key library " + path);
      }
      std::stringstream buf;
      buf << in.rdbuf();
      libs.push_back(topo::readLibrary(buf.str()));
    }
    const auto match = libs.size() == 1 ? topo::matchLibrary(fp, libs.front())
                                        : topo::matchLibraries(rows, libs, 7, colours);
    if (!perAtomPath.empty()) {
      // label index in sorted order of the labels seen, 0 for unmatched
      std::vector<std::string> labels;
      for (const auto &kv : match.counts) {
        if (!kv.first.empty()) {
          labels.push_back(kv.first);
        }
      }
      std::vector<int> values(static_cast<std::size_t>(cloud.nop), 0);
      for (std::size_t a = 0; a < match.labels.size() && a < values.size(); ++a) {
        const auto it = std::find(labels.begin(), labels.end(), match.labels[a]);
        values[a] = it == labels.end() ? 0 : static_cast<int>(it - labels.begin()) + 1;
      }
      writePerAtomDump(perAtomPath, cloud, "label", values);
      os << colorizer.longOption("labels");
      for (std::size_t i = 0; i < labels.size(); ++i) {
        os << " " << i + 1 << "=" << labels[i];
      }
      os << " ";
    }
    os << colorizer.heading("nop") << " " << rows.size() << " "
       << colorizer.longOption("graph") << " " << name << " "
       << colorizer.longOption("hops") << " " << hops << " "
       << colorizer.longOption("method") << " " << fp.method << " "
       << colorizer.longOption("matched") << " " << match.matched;
    for (const auto &kv : match.counts) {
      os << " " << (kv.first.empty() ? std::string("unmatched") : kv.first) << "=" << kv.second;
    }
    if (libs.size() > 1) {
      std::map<int, int> byDepth;
      for (const int d : match.depth) {
        if (d > 0) {
          byDepth[d] += 1;
        }
      }
      os << " " << colorizer.longOption("depth");
      for (auto it = byDepth.rbegin(); it != byDepth.rend(); ++it) {
        os << " " << it->first << "=" << it->second;
      }
    }
    os << "\n";
    return 0;
  }
  os << colorizer.heading("nop") << " " << rows.size() << " "
     << colorizer.longOption("graph") << " " << name << " "
     << colorizer.longOption("hops") << " " << hops << " "
     << colorizer.longOption("method") << " " << fp.method << " "
     << colorizer.longOption("key") << " " << fp.key << " "
     << colorizer.longOption("classes") << " " << fp.classes.size() << " "
     << colorizer.longOption("rings");
  for (std::size_t sz = 3; sz < fp.ringCensus.size(); ++sz) {
    os << " " << sz << ":" << fp.ringCensus[sz];
  }
  if (byAtom) {
    std::vector<int> species(rows.size(), 0);
    for (int i = 0; i < cloud.nop && static_cast<std::size_t>(i) < species.size(); ++i) {
      species[static_cast<std::size_t>(i)] = cloud.pts[static_cast<std::size_t>(i)].type;
    }
    const auto &rings = carry.rings(GraphSlot::single, rows, std::max(hops, 6));
    const auto table = cage::formerRows(rows, rings, species, -1);
    os << "\n";
    for (const auto &row : table) {
      os << "atom " << row.index << " species " << row.species << " coord " << row.coord
         << " homopolar " << row.homopolar << " rings";
      for (const auto &kv : row.rings) {
        os << " " << kv.first << ":" << kv.second;
      }
      os << "\n";
    }
    return 0;
  }
  // the most populated classes, largest first
  std::vector<std::pair<int, std::string>> top;
  for (const auto &kv : fp.classes) {
    top.emplace_back(kv.second, kv.first);
  }
  std::sort(top.begin(), top.end(), [](const auto &a, const auto &b) {
    return a.first != b.first ? a.first > b.first : a.second < b.second;
  });
  if (!perAtomPath.empty()) {
    // class rank by population, 1 the most populated
    std::map<std::string, int> rank;
    for (std::size_t t = 0; t < top.size(); ++t) {
      rank[top[t].second] = static_cast<int>(t) + 1;
    }
    std::vector<int> values(static_cast<std::size_t>(cloud.nop), 0);
    for (std::size_t a = 0; a < fp.atomKeys.size() && a < values.size(); ++a) {
      values[a] = rank[fp.atomKeys[a]];
    }
    writePerAtomDump(perAtomPath, cloud, "class", values);
  }
  os << " " << colorizer.longOption("top");
  for (std::size_t t = 0; t < top.size() && t < 5; ++t) {
    os << " " << top[t].second << "=" << top[t].first;
  }
  os << "\n";
  return 0;
}

// Ions read against the seeded cage assignment of the water: how many sit
// in ice, at the front, or in liquid by their first water shell.
int cmdIons(std::ostream &os, FrameCarry &carry, Cloud &cloud, double cutoff,
            int typeI, int k,
            bool complete, const std::vector<int> &ionTypes, double ionCutoff,
            const std::string &perAtomPath = {},
            const std::vector<int> &waterTypes = {}) {
  if (cloud.nop == 0) {
    os << colorizer.heading("nop") << " 0\n";
    return 0;
  }
  int typ = typeOf(cloud, typeI);
  if (!ionTypes.empty() &&
      std::find(ionTypes.begin(), ionTypes.end(), typ) != ionTypes.end()) {
    for (int i = 0; i < cloud.nop; ++i) {
      const int t = cloud.pts[static_cast<std::size_t>(i)].type;
      if (std::find(ionTypes.begin(), ionTypes.end(), t) == ionTypes.end()) {
        typ = t;
        break;
      }
    }
  }
  const double cand = cutoff + 1.5;
  auto graphs =
      waterTypes.empty()
          ? nneigh::kNearestNeighbourPair(cloud, k, cand, typ)
          : nneigh::kNearestNeighbourPair(cloud, k, cand, waterTypes);
  auto idxS = nneigh::neighbourListByIndex(cloud, graphs.first);
  auto idxU = nneigh::neighbourListByIndex(cloud, graphs.second);
  const auto sixS = sixRings(carry.rings(GraphSlot::strict, idxS, 6));
  const auto sixU = sixRings(carry.rings(GraphSlot::permissive, idxU, 6));
  const auto aff = ring::seededCageAffiliation(
      sixS, carry.affiliation(GraphSlot::strict, sixS, idxS), sixU,
      carry.affiliation(GraphSlot::permissive, sixU, idxU), idxU, complete);
  std::vector<bool> ice(static_cast<std::size_t>(cloud.nop), false);
  int nIce = 0;
  for (std::size_t i = 0; i < ice.size() && i < aff.hc.size(); ++i) {
    if (!isWaterType(cloud.pts[i].type, waterTypes)) {
      continue;
    }
    ice[i] = aff.hc[i] || aff.ddc[i];
    nIce += ice[i] ? 1 : 0;
  }
  std::vector<int> ions;
  for (int i = 0; i < cloud.nop; ++i) {
    const int t = cloud.pts[static_cast<std::size_t>(i)].type;
    if (std::find(ionTypes.begin(), ionTypes.end(), t) != ionTypes.end()) {
      ions.push_back(i);
    }
  }
  const int waterTyp = waterTypes.empty() ? typ : waterTypes.front();
  const auto env = site::ionEnvironment(cloud, ice, ions, waterTyp, ionCutoff);
  if (!perAtomPath.empty()) {
    // water: 0 liquid, 1 ice; ions: 2 liquid, 3 front, 4 ice
    std::vector<int> values(static_cast<std::size_t>(cloud.nop), 0);
    for (std::size_t i = 0; i < values.size() && i < ice.size(); ++i) {
      values[i] = ice[i] ? 1 : 0;
    }
    for (std::size_t i = 0; i < env.ion.size(); ++i) {
      values[static_cast<std::size_t>(env.ion[i])] = 2 + static_cast<int>(env.state[i]);
    }
    writePerAtomDump(perAtomPath, cloud, "state", values);
  }
  double meanShell = 0.0;
  double meanFraction = 0.0;
  for (std::size_t i = 0; i < env.ion.size(); ++i) {
    meanShell += env.shell[i];
    meanFraction += env.iceFraction[i];
  }
  if (!env.ion.empty()) {
    meanShell /= static_cast<double>(env.ion.size());
    meanFraction /= static_cast<double>(env.ion.size());
  }
  const auto census = site::iceClusterIonCensus(cloud, ice, idxU, ions, ionCutoff);
  int ionsInIce = 0;
  for (int n : census.ionsInCluster) {
    ionsInIce += n;
  }
  os << colorizer.heading("nop") << " " << cloud.nop << " "
     << colorizer.longOption("ice") << " " << nIce << " "
     << colorizer.longOption("ions") << " " << env.ion.size() << " "
     << colorizer.longOption("in-ice") << " " << env.nIce << " "
     << colorizer.longOption("front") << " " << env.nFront << " "
     << colorizer.longOption("liquid") << " " << env.nLiquid << " "
     << colorizer.longOption("shell") << " " << meanShell << " "
     << colorizer.longOption("shell-ice") << " " << meanFraction << " "
     << colorizer.longOption("clusters") << " " << census.nClusters << " "
     << colorizer.longOption("ions-in-ice") << " " << ionsInIce << "\n";
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

// CHILL and CHILL+ read one graph: cutoff (the default) or mutual. The
// cages default (seeded) falls back to cutoff.
std::vector<std::vector<int>>
chillGraph(Cloud &cloud, double cutoff, int typ, const std::string &graphName,
           const std::vector<int> &waterTypes) {
  const std::string name =
      (graphName.empty() || graphName == "seeded") ? std::string("cutoff") : graphName;
  const auto graph = nneigh::bondGraphFromName(name);
  if (graph == nneigh::BondGraph::Cutoff) {
    return nneigh::neighListO(cutoff, cloud,
                              waterTypes.empty() ? typ : waterTypes.front());
  }
  if (graph != nneigh::BondGraph::KnnMutual) {
    throw std::invalid_argument("CHILL and CHILL+ take --graph cutoff or mutual");
  }
  return waterGraph(cloud, cutoff, typ, waterTypes, 4, true);
}

int cmdChillPlus(std::ostream &os, FrameCarry &carry, Cloud &cloud,
                 double cutoff, int typeI,
                 const std::string &graphName, const std::string &ruleSpec,
                 bool layers = false,
                 bool tumLayers = false, int axis = 2, double layerWidth = 3.7,
                 const std::vector<int> &waterTypes = {}) {
  if (cloud.nop == 0) {
    printCounts(os, cloud);
    return 0;
  }
  const int typ = typeOf(cloud, typeI);
  auto nList = chillGraph(cloud, cutoff, typ, graphName, waterTypes);
  chill::applyChillPlus(cloud, nList, chill::chillPlusRuleFromFile(ruleSpec),
                        false);
  chill::clearIceTypesExcept(cloud, waterTypes);
  printCounts(os, cloud, waterTypes);
  if (layers) {
    const auto st = topoparam::layerCubicity(cloud, axis, layerWidth);
    os << colorizer.longOption("phi-c") << " " << st.phiC << " "
       << colorizer.longOption("stack") << " " << st.sequence << "\n";
  }
  if (tumLayers) {
    auto idx = nneigh::neighbourListByIndex(cloud, nList);
    const auto six = sixRings(carry.rings(GraphSlot::single, idx, 6));
    const auto planes = ring::stackingPlanes(six, idx);
    const auto st =
        topoparam::tumLayerStack(cloud, six, planes.basal, planes.equatorial,
                                 axis, layerWidth);
    os << colorizer.longOption("tum-phi-c") << " " << st.phiC << " "
       << colorizer.longOption("tum-stack") << " " << st.sequence << "\n";
  }
  return 0;
}

int cmdChill(std::ostream &os, Cloud &cloud, double cutoff, int typeI,
             const std::string &graphName,
             const std::vector<int> &waterTypes = {}) {
  if (cloud.nop == 0) {
    printCounts(os, cloud);
    return 0;
  }
  const int typ = typeOf(cloud, typeI);
  auto nList = chillGraph(cloud, cutoff, typ, graphName, waterTypes);
  chill::getCorrel(cloud, nList, false);
  chill::getIceTypeNoPrint(cloud, nList, false);
  chill::clearIceTypesExcept(cloud, waterTypes);
  printCounts(os, cloud, waterTypes);
  return 0;
}

int cmdF4(std::ostream &os, Cloud &cloud, double cutoff, int typeI, int hType) {
  if (cloud.nop == 0) {
    os << colorizer.heading("nop") << " 0 "
       << colorizer.longOption("f4") << " nan\n";
    return 0;
  }
  int typ = typeI;
  if (typ <= 0) {
    for (const auto &p : cloud.pts) {
      if (p.type != hType) {
        typ = p.type;
        break;
      }
    }
  }
  if (typ <= 0) {
    typ = typeOf(cloud, typeI);
  }
  const double cand = cutoff + 1.5;
  auto nList = nneigh::kNearestNeighbourList(cloud, 4, cand, typ, true);
  const auto f4 = topoparam::rodgerF4(cloud, nList, typ, hType);
  int nFinite = 0;
  for (double v : f4) {
    nFinite += std::isfinite(v) ? 1 : 0;
  }
  os << colorizer.heading("nop") << " " << cloud.nop << " "
     << colorizer.longOption("f4") << " " << topoparam::meanFinite(f4) << " "
     << colorizer.longOption("n") << " " << nFinite << "\n";
  return 0;
}

int cmdSteinhardt(std::ostream &os, Cloud &cloud, double cutoff, int typeI) {
  if (cloud.nop == 0) {
    os << colorizer.heading("nop") << " 0\n";
    return 0;
  }
  const int typ = typeOf(cloud, typeI);
  const auto q3 = chill::steinhardtQlCutoff(cloud, cutoff, 3, typ, true);
  const auto q12 = chill::steinhardtQlCutoff(cloud, cutoff, 12, typ, true);
  os << colorizer.heading("nop") << " " << cloud.nop << " "
     << colorizer.longOption("q3") << " " << topoparam::meanFinite(q3.ql) << " "
     << colorizer.longOption("q3bar") << " " << topoparam::meanFinite(q3.qlBar)
     << " " << colorizer.longOption("q12") << " "
     << topoparam::meanFinite(q12.ql) << " "
     << colorizer.longOption("q12bar") << " "
     << topoparam::meanFinite(q12.qlBar) << "\n";
  return 0;
}

// Several signatures joined with '+' share one ring network and one search.
std::vector<cage::Signature> parseSignatures(const std::string &spec) {
  std::vector<cage::Signature> sigs;
  std::string part;
  std::istringstream in(spec);
  while (std::getline(in, part, '+')) {
    sigs.push_back(cage::Signature::parse(part));
  }
  return sigs;
}

int cmdCages(std::ostream &os, FrameCarry &carry, Cloud &cloud, double cutoff,
             int typeI, int k,
             const std::string &graphName, bool complete = false,
             const std::string &signatureSpec = {},
             const std::vector<int> &guestTypes = {}, double guestRadius = 4.0,
             const std::string &perAtomPath = {}, bool inside = false,
             bool incomplete = false, const std::vector<int> &waterTypes = {},
             bool shape = false, const std::string &cageEdge = "2.65,3.15") {
  if (cloud.nop == 0) {
    if (!signatureSpec.empty()) {
      for (const auto &sig : parseSignatures(signatureSpec)) {
        os << colorizer.heading("nop") << " 0 "
           << colorizer.longOption("graph") << " " << graphName << " "
           << colorizer.longOption("signature") << " " << sig.str()
           << " cages 0 atoms 0\n";
      }
      return 0;
    }
    os << colorizer.heading("nop") << " 0 "
       << colorizer.longOption("graph") << " " << graphName << " "
       << iceColor("hexagonal") << " 0 " << iceColor("cubic") << " 0 "
       << iceColor("both") << " 0 " << iceColor("water") << " 0\n";
    return 0;
  }
  const int typ = typeOf(cloud, typeI);
  const double cand = cutoff + 1.5;
  auto knnOf = [&](bool mutual) {
    if (!waterTypes.empty()) {
      return nneigh::kNearestNeighbourList(cloud, k, cand, waterTypes, mutual);
    }
    return nneigh::kNearestNeighbourList(cloud, k, cand, typ, mutual);
  };
  auto knnPairOf = [&]() {
    if (!waterTypes.empty()) {
      return nneigh::kNearestNeighbourPair(cloud, k, cand, waterTypes);
    }
    return nneigh::kNearestNeighbourPair(cloud, k, cand, typ);
  };

  if (!signatureSpec.empty()) {
    const auto sigs = parseSignatures(signatureSpec);
    const auto emit = [&](const std::vector<std::vector<int>> &idx) {
      int depth = 3;
      for (const auto &sig : sigs) {
        depth = std::max(depth, sig.maxRingSize());
      }
      const auto rings = ring::cageRings(
          carry.rings(GraphSlot::single, idx, depth), cloud, cutoff);
      const auto rank = ring::positionRank(cloud);
      const auto window = cage::parseEdgeWindow(cageEdge);
      std::vector<std::array<double, 3>> xyz;
      const std::vector<std::array<double, 3>> *xyzPtr = nullptr;
      const std::vector<double> *boxPtr = nullptr;
      if (window.enabled) {
        xyz.resize(static_cast<std::size_t>(cloud.nop));
        for (int i = 0; i < cloud.nop; ++i) {
          const auto &p = cloud.pts[static_cast<std::size_t>(i)];
          xyz[static_cast<std::size_t>(i)] = {p.x, p.y, p.z};
        }
        xyzPtr = &xyz;
        boxPtr = &cloud.box;
      }
      const auto foundOf =
          cage::findBySignatures(rings, idx, sigs, rank, xyzPtr, boxPtr, window);
      for (size_t si = 0; si < sigs.size(); ++si) {
        const auto &sig = sigs[si];
        const auto &found = foundOf[si];
        std::vector<cage::FoundCage> cups;
        if (incomplete) {
          cups = cage::findIncompleteBySignature(rings, sig, 0);
        }
        std::set<int> atoms;
        for (const auto &c : found) {
          atoms.insert(c.vertices.begin(), c.vertices.end());
        }
        os << colorizer.heading("nop") << " " << cloud.nop << " "
           << colorizer.longOption("graph") << " " << graphName << " "
           << colorizer.longOption("signature") << " " << sig.str() << " cages "
           << found.size() << " atoms " << atoms.size();
        if (incomplete) {
          os << " " << colorizer.longOption("incomplete") << " " << cups.size();
        }
        if (!guestTypes.empty()) {
          // guests (methane, THF, ions) placed in the found cages by the
          // periodic centroid of each cage's vertices
          std::vector<int> guests;
          for (int i = 0; i < cloud.nop; ++i) {
            const int t = cloud.pts[static_cast<std::size_t>(i)].type;
            if (std::find(guestTypes.begin(), guestTypes.end(), t) != guestTypes.end()) {
              guests.push_back(i);
            }
          }
          std::vector<std::vector<int>> cages;
          cages.reserve(found.size());
          for (const auto &c : found) {
            cages.push_back(c.vertices);
          }
          site::GuestOccupancy occ;
          if (inside) {
            std::vector<std::vector<int>> faces;
            faces.reserve(found.size());
            for (const auto &c : found) {
              faces.push_back(c.faces);
            }
            occ = site::guestOccupancyInside(cloud, rings, faces, guests);
          } else {
            occ = site::guestOccupancy(cloud, cages, guests, guestRadius);
          }
          os << " " << colorizer.longOption("guests") << " " << guests.size() << " occupied "
             << occ.occupied << " multiple " << occ.multiply << " free " << occ.free
             << (inside ? " inside 1" : "");
          if (!occ.occupancyHistogram.empty()) {
            os << " " << colorizer.longOption("occ");
            for (std::size_t k = 0; k < occ.occupancyHistogram.size(); ++k) {
              os << " " << k << ":" << occ.occupancyHistogram[k];
            }
          }
        }
        if (shape) {
          std::vector<std::array<double, 3>> all(static_cast<std::size_t>(cloud.nop));
          for (int i = 0; i < cloud.nop; ++i) {
            const auto &p = cloud.pts[static_cast<std::size_t>(i)];
            all[static_cast<std::size_t>(i)] = {p.x, p.y, p.z};
          }
          os << " " << colorizer.longOption("shape");
          for (const auto &c : found) {
            const auto xyz = cage::coordsOfVertices(all, c.vertices);
            const auto sh = cage::shapeOfVertices(xyz);
            os << " " << sh.nVertices << ":"
               << (sh.status == 0 ? sh.pointGroup : "off");
          }
        }
        os << "\n";
      }
      return 0;
    };
    if (graphName == "seeded") {
      auto graphs = knnPairOf();
      auto idxU = nneigh::neighbourListByIndex(cloud, graphs.second);
      return emit(idxU);
    }
    const auto graph = nneigh::bondGraphFromName(graphName);
    std::vector<std::vector<int>> nList;
    if (graph == nneigh::BondGraph::Cutoff) {
      nList = nneigh::neighListO(cutoff, cloud,
                                waterTypes.empty() ? typ : waterTypes.front());
    } else {
      const bool mutual = graph == nneigh::BondGraph::KnnMutual;
      nList = knnOf(mutual);
    }
    auto idx = nneigh::neighbourListByIndex(cloud, nList);
    return emit(idx);
  }

  int ih = 0;
  int ic = 0;
  int both = 0;
  int water = 0;
  const auto tallyAtoms = [&](const std::vector<bool> &hc,
                              const std::vector<bool> &ddc) {
    // Disjoint: hexagonal only, cubic only, both, water. A molecule on
    // both cage types is `both`, not hexagonal. Per-atom codes are
    // 0 water, 1 hexagonal, 2 cubic, 3 both. The four counters sum to
    // the number of water molecules.
    ih = ic = both = water = 0;
    const int n = static_cast<int>(hc.size());
    for (int i = 0; i < n; ++i) {
      if (i < cloud.nop && !isWaterType(cloud.pts[static_cast<std::size_t>(i)].type,
                                        waterTypes)) {
        continue;
      }
      const bool h = hc[static_cast<std::size_t>(i)];
      const bool d = ddc[static_cast<std::size_t>(i)];
      if (h && d) {
        ++both;
      } else if (h) {
        ++ih;
      } else if (d) {
        ++ic;
      } else {
        ++water;
      }
    }
    if (!perAtomPath.empty()) {
      // 0 water, 1 hexagonal cage, 2 double-diamond cage, 3 both
      std::vector<int> values(static_cast<std::size_t>(cloud.nop), 0);
      for (int i = 0; i < n && i < cloud.nop; ++i) {
        const bool h = hc[static_cast<std::size_t>(i)];
        const bool d = ddc[static_cast<std::size_t>(i)];
        values[static_cast<std::size_t>(i)] = h && d ? 3 : h ? 1 : d ? 2 : 0;
      }
      writePerAtomDump(perAtomPath, cloud, "cage", values);
    }
  };

  if (graphName == "seeded") {
    auto graphs = knnPairOf();
    const auto &mutual = graphs.first;
    const auto &uni = graphs.second;
    auto idxS = nneigh::neighbourListByIndex(cloud, mutual);
    auto idxU = nneigh::neighbourListByIndex(cloud, uni);
    // tum::cageCounts reads one graph and has no ring completion, so the
    // seeded assignment stays on the host whatever SEAMS_OFFLOAD says.
    // Hash first. A match is checked again by comparing every edge.
    if (seams::provenEqual(seams::edgeSetOf(idxS), seams::edgeSetOf(idxU))) {
      const auto six = ring::contractibleRings(
          sixRings(carry.rings(GraphSlot::strict, idxS, 6)), cloud);
      const auto &once = carry.affiliation(GraphSlot::strict, six, idxS);
      const auto aff =
          ring::seededCageAffiliation(six, once, six, once, idxS, complete);
      tallyAtoms(aff.hc, aff.ddc);
    } else {
      const auto sixS = ring::contractibleRings(
          sixRings(carry.rings(GraphSlot::strict, idxS, 6)), cloud);
      const auto sixU = ring::contractibleRings(
          sixRings(carry.rings(GraphSlot::permissive, idxU, 6)), cloud);
      const auto aff = ring::seededCageAffiliation(
          sixS, carry.affiliation(GraphSlot::strict, sixS, idxS), sixU,
          carry.affiliation(GraphSlot::permissive, sixU, idxU), idxU,
          complete);
      tallyAtoms(aff.hc, aff.ddc);
    }
  } else {
    const auto graph = nneigh::bondGraphFromName(graphName);
    std::vector<std::vector<int>> nList;
    if (graph == nneigh::BondGraph::Cutoff) {
      nList = nneigh::neighListO(cutoff, cloud,
                                waterTypes.empty() ? typ : waterTypes.front());
    } else {
      const bool mutual = graph == nneigh::BondGraph::KnnMutual;
      nList = knnOf(mutual);
    }
    auto idx = nneigh::neighbourListByIndex(cloud, nList);
    if (tum::preferOffload()) {
      const auto counts = tum::cageCounts(idx);
      std::vector<bool> hc(static_cast<std::size_t>(cloud.nop), false);
      std::vector<bool> ddc(static_cast<std::size_t>(cloud.nop), false);
      const int n = std::min(cloud.nop, static_cast<int>(counts.atomHc.size()));
      for (int i = 0; i < n; ++i) {
        hc[static_cast<std::size_t>(i)] = counts.atomHc[static_cast<std::size_t>(i)] != 0;
        ddc[static_cast<std::size_t>(i)] = counts.atomDdc[static_cast<std::size_t>(i)] != 0;
      }
      tallyAtoms(hc, ddc);
    } else {
      const auto six = ring::contractibleRings(
          sixRings(carry.rings(GraphSlot::single, idx, 6)), cloud);
      const auto &aff = carry.affiliation(GraphSlot::single, six, idx);
      std::vector<bool> hc(static_cast<std::size_t>(cloud.nop), false);
      std::vector<bool> ddc(static_cast<std::size_t>(cloud.nop), false);
      for (std::size_t r = 0; r < six.size(); ++r) {
        for (const int a : six[r]) {
          hc[static_cast<std::size_t>(a)] = hc[static_cast<std::size_t>(a)] || aff.hc[r];
          ddc[static_cast<std::size_t>(a)] = ddc[static_cast<std::size_t>(a)] || aff.ddc[r];
        }
      }
      tallyAtoms(hc, ddc);
    }
  }
  os << colorizer.heading("nop") << " " << waterNop(cloud, waterTypes) << " "
     << colorizer.longOption("graph") << " " << graphName << " "
     << iceColor("hexagonal") << " " << ih << " " << iceColor("cubic") << " "
     << ic << " " << iceColor("both") << " " << both << " "
     << iceColor("water") << " " << water << "\n";
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

bool parseMolList(std::string_view value, std::vector<int> &ids) {
  ids.clear();
  if (value.empty()) {
    return false;
  }
  std::size_t start = 0;
  while (start <= value.size()) {
    const std::size_t comma = value.find(',', start);
    const auto token =
        trimView(value.substr(start, comma == std::string_view::npos
                                         ? std::string_view::npos
                                         : comma - start));
    if (token.empty()) {
      return false;
    }
    const auto dash = token.find('-');
    try {
      if (dash == std::string_view::npos) {
        ids.push_back(parseIntegral<int>(token));
      } else {
        const int lo = parseIntegral<int>(token.substr(0, dash));
        const int hi = parseIntegral<int>(token.substr(dash + 1));
        if (hi < lo) {
          return false;
        }
        for (int mol = lo; mol <= hi; ++mol) {
          ids.push_back(mol);
        }
      }
    } catch (const ParsingException &) {
      return false;
    }
    if (comma == std::string_view::npos) {
      break;
    }
    start = comma + 1;
  }
  return !ids.empty();
}

bool scoringCommand(std::string_view cmd) {
  return cmd == "chill" || cmd == "chill-plus" || cmd == "chill_plus" ||
         cmd == "cages" || cmd == "fingerprint" || cmd == "steinhardt";
}

} // namespace

int main(int argc, char *argv[]) {
  colorizer = colorizerForFile(environmentColorStatus(), stdout);

  auto cfg = seams::cfg::load(seams::cfg::pathFromArgv(argc, argv));
  int frame = cfg.frame;
  int last = cfg.last;
  int jobs = cfg.jobs;
  bool carryFrames = cfg.carry;
  int typeI = cfg.type;
  double cutoff = cfg.cutoff;
  std::string cageEdge = cfg.cageEdge;
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
  bool completeFlag = false;
  std::string signatureSpec;
  std::string guestTypesFlag;
  double guestRadius = 4.0;
  std::string perAtomPath;
  bool colourTypes = false;
  std::string libraryPath;
  std::string emitLabel;
  int hops = 2;
  std::string ionTypesFlag;
  double ionCutoff = 0.0;
  bool insideFlag = false;
  bool layersFlag = false;
  bool tumLayersFlag = false;
  bool incompleteFlag = false;
  bool shapeFlag = false;
  bool byAtomFlag = false;
  std::string waterTypesFlag;
  std::string molFlag;
  bool molWater = false;
  bool molSet = false;
  std::vector<int> molIds;
  std::string subsetFlag;
  int rdfTypeI = 0;
  int rdfTypeJ = 0;
  bool typesSet = false;
  int bins = 0;
  int densAxis = 2;
  std::string axisFlag;
  std::string outputFormat = "text";
  std::string chillRule;
  bool chillRuleSet = false;

  const char *progname = argv[0];
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
                 .help("Last frame (inclusive). Omit for a single --frame"
#ifdef SEAMS_HAS_MPI
                       ". Under mpiexec the ranks share a dump's frames, "
                       "each with --jobs workers; SEAMS_MPI=0 keeps MPI off"
#endif
                       )
                 .handler([&](std::string_view value) {
                   last = parseIntegral<int>(value);
                 }));

  parser.add(Option("--jobs", "-j")
                 .argName("N")
                 .help("Parallel frame workers (OpenMP). 1 is serial")
                 .handler([&](std::string_view value) {
                   jobs = parseIntegral<int>(value);
                 }));

  parser.add(Option("--rebuild")
                 .help("Recompute rings and cages on every frame. A one-worker "
                       "frame range otherwise updates them from the previous "
                       "frame; the output is the same")
                 .handler([&]() { carryFrames = false; }));

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

  parser.add(Option("--cage-edge")
                 .argName("SPEC")
                 .help("cages --signature: keep a face only when every edge "
                       "lies in LO,HI (coordinate units of the file; default "
                       "2.65,3.15). none uses the face census alone "
                       "(SEAMS_CAGE_EDGE)")
                 .handler([&](std::string_view value) {
                   cageEdge = std::string(value);
                 }));

  parser.add(Option("--pair-types")
                 .argName("I,J")
                 .help("Pair types for site-site rdf/cn (default I=J=--type)")
                 .handler([&](std::string_view value) {
                   typesFlag = std::string(value);
                 }));

  // Exact alias. Argum returns an exact key before any longer name
  // that starts with the same characters, so --type and --types both
  // bind. A shorter token such as --ty matches both and is ambiguous.
  parser.add(Option("--types")
                 .argName("I,J")
                 .help("Deprecated exact alias of --pair-types")
                 .handler([&](std::string_view value) {
                   std::cerr << colorizer.warning(
                                    "--types is deprecated; use --pair-types")
                             << "\n";
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

  parser.add(Option("--rule")
                 .argName("NAME")
                 .help("chill-plus rule: chill+ (published default) or "
                       "ovito-chill+. A path to a readable file is scanned "
                       "for a rule or chill-plus key. Comments are skipped. "
                       "The file is not executed as Lua or YAML")
                 .handler([&](std::string_view value) {
                   chillRule = std::string(value);
                   chillRuleSet = true;
                 }));

  parser.add(Option("--graph")
                 .argName("KIND")
                 .help("Bond graph for cages: cutoff | knn | knn-union | seeded; "
                       "chill, chill-plus: cutoff (default) | mutual")
                 .handler([&](std::string_view value) { graph = value; }));

  parser.add(Option("--complete")
                 .help("Seeded cages: fill the last vertex of six-rings whose "
                       "other vertices carry a label (ring completion)")
                 .handler([&]() { completeFlag = true; }));

  parser.add(Option("--shape")
                 .help("cages --signature: point group of each found cage, "
                       "from those vertices only")
                 .handler([&]() { shapeFlag = true; }));

  parser.add(Option("--signature")
                 .argName("SPEC")
                 .help("cages: ring-size census (4:6,6:8) or a named table "
                       "entry (sodalite|alpha|512|51262|51264|51268|sh|hc|ddc); "
                       "join several with + for one line each")
                 .handler([&](std::string_view value) {
                   signatureSpec = std::string(value);
                 }));

  parser.add(Option("--guest-types")
                 .argName("T,U")
                 .help("cages --signature: LAMMPS types of guests (methane, THF, ions) "
                       "placed in the found cages by their periodic centroids")
                 .handler([&](std::string_view value) {
                   guestTypesFlag = std::string(value);
                 }));

  parser.add(Option("--per-atom")
                 .argName("FILE")
                 .help("cages|fingerprint|ions: append a LAMMPS dump frame with one extra "
                       "per-atom column (cage, class or label, state) for OVITO or VMD")
                 .handler([&](std::string_view value) {
                   perAtomPath = std::string(value);
                 }));

  parser.add(Option("--guest-radius")
                 .argName("R")
                 .help("cages --signature: a guest belongs to the nearest cage centre "
                       "within R (default 4.0)")
                 .handler([&](std::string_view value) {
                   guestRadius = parseFloatingPoint<double>(value);
                 }));

  parser.add(Option("--inside")
                 .help("cages --guest-types: assign by ray-parity inside the "
                       "fan-triangulated faces, not nearest centroid")
                 .handler([&]() { insideFlag = true; }));

  parser.add(Option("--layers")
                 .help("chill-plus: literature I_sd reference from CHILL+ "
                       "molecule bins (Phi_c and H/C string)")
                 .handler([&]() { layersFlag = true; }));

  parser.add(Option("--tum-layers")
                 .help("chill-plus: TUM stacking from HC-basal and "
                       "DDC-equatorial rings, not CHILL+ molecules")
                 .handler([&]() { tumLayersFlag = true; }));

  parser.add(Option("--water-types")
                 .argName("T,U")
                 .help("Neighbour graph on these LAMMPS types only (OW). "
                       "Substrate and ions never enter the 4-NN list. "
                       "A type that is not in the frame is an error")
                 .handler([&](std::string_view value) {
                   waterTypesFlag = std::string(value);
                 }));

  parser.add(Option("--mol")
                 .argName("SPEC")
                 .help("Water oxygens by molecule: water, or ids such as "
                       "1,5-12. Neighbours, CHILL, CHILL+, rings and cages "
                       "then ignore ions. A missing id is an error")
                 .handler([&](std::string_view value) {
                   molFlag = std::string(value);
                 }));

  parser.add(Option("--incomplete")
                 .help("cages --signature: also report cups and incomplete "
                       "cages on the ring graph")
                 .handler([&]() { incompleteFlag = true; }));

  parser.add(Option("--colour-types")
                 .help("fingerprint: colour vertices by LAMMPS type, so species "
                       "never match across types")
                 .handler([&]() { colourTypes = true; }));

  parser.add(Option("--library")
                 .argName("FILE")
                 .help("fingerprint: key library (from --emit-library) to name atoms "
                       "by; several, comma separated, at different hop counts name "
                       "each atom by the deepest that knows it")
                 .handler([&](std::string_view value) {
                   libraryPath = std::string(value);
                 }));

  parser.add(Option("--emit-library")
                 .argName("LABEL")
                 .help("fingerprint: print the frame's distinct keys as library lines "
                       "under LABEL instead of the summary")
                 .handler([&](std::string_view value) {
                   emitLabel = std::string(value);
                 }));

  parser.add(Option("--hops")
                 .argName("N")
                 .help("fingerprint: bonds from the centre in each local key (default 2)")
                 .handler([&](std::string_view value) {
                   hops = parseIntegral<int>(value);
                 }));

  parser.add(Option("--by-atom")
                 .help("fingerprint: one row per atom with coordination, "
                       "same-species bonds, and primitive rings through that atom")
                 .handler([&]() { byAtomFlag = true; }));

  parser.add(Option("--ion-types")
                 .argName("I,J")
                 .help("ions: LAMMPS types read against the cage assignment")
                 .handler([&](std::string_view value) {
                   ionTypesFlag = std::string(value);
                 }));

  parser.add(Option("--ion-cutoff")
                 .argName("ANGSTROM")
                 .help("ions: first water shell radius (default --cutoff)")
                 .handler([&](std::string_view value) {
                   ionCutoff = parseFloatingPoint<double>(value);
                 }));

  parser.add(Option("--htype")
                 .argName("I")
                 .help("Hydrogen atom type for hbonds and f4")
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
                 .help("read | chill | chill-plus | ovito-chill+ | cages | symmetry | "
                       "fingerprint | ions | "
                       "f4 | steinhardt | rdf | cn | hbonds | pairs | density-z | "
                       "domains")
                 .occurs(zeroOrOneTime)
                 .handler([&](std::string_view value) { cmd = value; }));

  parser.add(Positional("file")
                 .help("Trajectory file")
                 .occurs(zeroOrOneTime)
                 .handler([&](std::string_view value) { file = value; }));

  try {
    parser.parse(argc, argv);
  } catch (...) {
    std::cerr << colorizer.error(exceptionText()) << "\n";
    std::cerr << colorizer.warning(parser.formatUsage(progname, colorizer))
              << "\n";
    return 2;
  }
#ifdef SEAMS_HAS_MPI
  // After parsing: --help, --version and --features leave through std::exit
  const FrameRanks ranks(&argc, &argv);
  if (ranks.active && environmentColorStatus() != ColorStatus::required) {
    // Open MPI gives ranks a terminal even when mpiexec's output is a file
    colorizer = {};
  }
  if (!ranks.sameCommand || (ranks.size > 1 && !perAtomPath.empty())) {
    if (ranks.rank == 0) {
      std::cerr << colorizer.error(
                       ranks.sameCommand
                           ? "--per-atom writes one file; run it on one rank"
                           : "every rank of a launch must run the same command "
                             "in the same place; SEAMS_MPI=0 runs each alone")
                << "\n";
    }
    return 2;
  }
#endif

  if (outputFormat != "text" && outputFormat != "json") {
    std::cerr << colorizer.error("bad --format (want text|json)") << "\n";
    return 2;
  }

  cfg.frame = frame;
  cfg.last = last;
  cfg.jobs = jobs;
  cfg.carry = carryFrames;
  cfg.type = typeI;
  cfg.cutoff = cutoff;
  try {
    cfg.cageEdge = cage::edgeWindowSpec(cage::parseEdgeWindow(cageEdge));
  } catch (...) {
    std::cerr << colorizer.error(exceptionText()) << "\n";
    std::cerr << colorizer.warning(parser.formatUsage(progname, colorizer))
              << "\n";
    return 2;
  }
  cfg.k = k;
  cfg.graph = graph;
  if (familyFlag) {
    try {
      family = site::parseFamily(familyNameFlag);
    } catch (...) {
      std::cerr << colorizer.error(exceptionText()) << "\n";
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
    try {
      if (!parseTypePair(typesFlag, rdfTypeI, rdfTypeJ)) {
        std::cerr << colorizer.error("bad --pair-types (want I,J)") << "\n";
        std::cerr << colorizer.warning(parser.formatUsage(progname, colorizer))
                  << "\n";
        return 2;
      }
    } catch (...) {
      std::cerr << colorizer.error(exceptionText()) << "\n";
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

  if (cmd == "ovito-chill+" || cmd == "ovito_chill_plus") {
    if (!chillRuleSet) {
      chillRule = "ovito-chill+";
    }
    cmd = "chill-plus";
  }
  if (!chillRuleSet && chillRule.empty()) {
    if (const char *fromEnv = std::getenv("CHILL_PLUS_RULE")) {
      chillRule = fromEnv;
    }
  }
  const bool iceCmd = cmd == "chill" || cmd == "chill-plus" ||
                      cmd == "chill_plus" || cmd == "cages" || cmd == "ions";
  if (iceCmd && !site::iceScoreAllowed(family)) {
    std::cerr << colorizer.error(site::refuseIceScore(family)) << "\n";
    return 2;
  }

  site::Table siteTable;
  siteTable.family = family;
  if (!siteSpec.empty()) {
    try {
      siteTable = parseSiteSpec(siteSpec, family);
    } catch (...) {
      std::cerr << colorizer.error(exceptionText()) << "\n";
      return 2;
    }
  }
  if (cmd == "ions" && ionTypesFlag.empty()) {
    std::cerr << colorizer.error("ions needs --ion-types") << "\n";
    return 2;
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
  std::vector<int> waterTypes;
  {
    std::string tok;
    std::istringstream in(waterTypesFlag);
    while (std::getline(in, tok, ',')) {
      if (!tok.empty()) {
        waterTypes.push_back(parseIntegral<int>(tok));
      }
    }
  }
  if (!molFlag.empty()) {
    const auto spec = trimView(molFlag);
    if (spec == "water") {
      molWater = true;
      molSet = true;
    } else if (!parseMolList(spec, molIds)) {
      std::cerr << colorizer.error(
                       "bad --mol (want water or ids such as 1,5-12)")
                << "\n";
      return 2;
    } else {
      molSet = true;
    }
  }
  std::vector<int> ionTypes;
  {
    std::string tok;
    std::istringstream in(ionTypesFlag);
    while (std::getline(in, tok, ',')) {
      if (!tok.empty()) {
        try {
          ionTypes.push_back(parseIntegral<int>(tok));
        } catch (const ParsingException &) {
          break;
        }
      }
    }
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

  auto runOne = [&](std::ostream &os, FrameCarry &carry, Cloud &cloud) {
    if (strictInput && cloud.nop == 0) {
      os << "input frame " << cloud.currentFrame << " contains no atoms\n";
      return 2;
    }
    const bool scoring = scoringCommand(cmd);
    std::vector<water::IonSite> brineIons;
    if (scoring && cloud.nop > 0 &&
        (molSet || !waterTypes.empty() ||
         (typeI <= 0 && frameSource(file) == FrameSource::dump))) {
      try {
        const bool guests = cmd == "cages" && !guestTypesFlag.empty();
        if (!waterTypes.empty() && !molSet) {
          water::Request request;
          request.types = waterTypes;
          request.ion_types = ionTypes;
          // Validates the types. The existing graph path still decides
          // which of those types bond. A missing type is the error.
          (void)water::selectOxygens(cloud, request);
        } else if (molSet ||
                   (typeI <= 0 && water::hasMultiAtomMolecule(cloud))) {
          water::Request request;
          request.ion_types = ionTypes;
          if (!waterTypes.empty()) {
            request.types = waterTypes;
          } else if (typeI > 0) {
            request.type = typeI;
          }
          if (molWater || (!molSet && typeI <= 0)) {
            request.mol_mode = water::MolMode::Water;
          } else if (molSet) {
            request.mol_mode = water::MolMode::Ids;
            request.molecules = molIds;
          }
          const auto selection = water::selectOxygens(cloud, request);
          brineIons = water::collectIons(cloud, selection);
          if (guests) {
            if (waterTypes.empty()) {
              waterTypes = selection.oxygen_types;
            }
          } else {
            cloud = water::subsetCloud(cloud, selection);
            waterTypes.clear();
          }
        } else if (typeI <= 0 && frameSource(file) == FrameSource::dump) {
          cloud = water::keepHistoricalOxygen(cloud);
        }
      } catch (...) {
        os << colorizer.error(exceptionText()) << "\n";
        return 2;
      }
    }
    // --family sets the flag in the same handler that replaces the
    // default, so an unset flag is waterIce.
    if (!familyFlag && uniqueTypeCount(cloud) > 2) {
      std::cerr << colorizer.warning(
                       "hint: more than two LAMMPS types; set --family "
                       "(default waterIce)")
                << "\n";
    }
    if (cmd == "read") {
      return cmdRead(os, cloud);
    }
    if (cmd == "symmetry") {
      try {
        return cmdSymmetry(os, cloud);
      } catch (const std::exception &e) {
        os << colorizer.error(e.what()) << "\n";
        return 2;
      }
    }
    if (cmd == "chill-plus" || cmd == "chill_plus") {
      try {
        const int rc =
            cmdChillPlus(os, carry, cloud, cutoff, typeI, graph, chillRule,
                         layersFlag, tumLayersFlag, densAxis, 3.7, waterTypes);
        if (rc == 0 && !brineIons.empty()) {
          const auto report = water::ionFront(cloud, brineIons);
          os << "ion-ice " << report.ice << " ion-liquid " << report.liquid
             << " ion-front " << report.front << "\n";
        }
        return rc;
      } catch (...) {
        os << colorizer.error(exceptionText()) << "\n";
        return 2;
      }
    }
    if (cmd == "chill") {
      try {
        const int rc = cmdChill(os, cloud, cutoff, typeI, graph, waterTypes);
        if (rc == 0 && !brineIons.empty()) {
          const auto report = water::ionFront(cloud, brineIons);
          os << "ion-ice " << report.ice << " ion-liquid " << report.liquid
             << " ion-front " << report.front << "\n";
        }
        return rc;
      } catch (...) {
        os << colorizer.error(exceptionText()) << "\n";
        return 2;
      }
    }
    if (cmd == "cages") {
      try {
        std::vector<int> guestTypes;
        std::string tok;
        std::istringstream in(guestTypesFlag);
        while (std::getline(in, tok, ',')) {
          if (!tok.empty()) {
            guestTypes.push_back(parseIntegral<int>(tok));
          }
        }
        return cmdCages(os, carry, cloud, cutoff, typeI, k, graph, completeFlag,
                        signatureSpec, guestTypes, guestRadius, perAtomPath,
                        insideFlag, incompleteFlag, waterTypes, shapeFlag,
                        cfg.cageEdge);
      } catch (...) {
        os << colorizer.error(exceptionText()) << "\n";
        return 2;
      }
    }
    if (cmd == "fingerprint") {
      try {
        return cmdFingerprint(os, carry, cloud, cutoff, typeI, k, graph, hops,
                              colourTypes,
                              libraryPath, emitLabel, perAtomPath, byAtomFlag);
      } catch (...) {
        os << colorizer.error(exceptionText()) << "\n";
        return 2;
      }
    }
    if (cmd == "ions") {
      std::vector<int> ionTypes;
      std::string tok;
      std::istringstream in(ionTypesFlag);
      while (std::getline(in, tok, ',')) {
        if (!tok.empty()) {
          ionTypes.push_back(parseIntegral<int>(tok));
        }
      }
      return cmdIons(os, carry, cloud, cutoff, typeI, k, completeFlag, ionTypes,
                     ionCutoff > 0.0 ? ionCutoff : cutoff, perAtomPath,
                     waterTypes);
    }
    if (cmd == "f4") {
      return cmdF4(os, cloud, cutoff, typeI, htype);
    }
    if (cmd == "steinhardt") {
      return cmdSteinhardt(os, cloud, cutoff, typeI);
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

  const bool scoreGuess =
      scoringCommand(cmd) && typeI <= 0 && waterTypes.empty() && !molSet &&
      frameSource(file) == FrameSource::dump;
  const bool loadAll =
      cmd == "pairs" || cmd == "ions" || cmd == "f4" ||
      (cmd == "cn" && ionsFlag) ||
      ((cmd == "rdf" || cmd == "cn") && typesSet && rdfTypeI != rdfTypeJ) ||
      (cmd == "density-z" && typeI <= 0) || cmd == "domains" ||
      !waterTypes.empty() || molSet || scoreGuess;
  const int loadType = loadAll ? -1 : typeI;

  if (last <= 0 || last == frame) {
#ifdef SEAMS_HAS_MPI
    if (ranks.rank != 0) {
      return 0;
    }
#endif
    Cloud cloud;
    try {
      cloud = load(file, frame, loadType);
    } catch (...) {
      std::cerr << colorizer.error(exceptionText()) << "\n";
      return 2;
    }
    std::ostringstream line;
    FrameCarry once(false);
    const int rc = runOne(line, once, cloud);
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
  // Set below, before the walk starts, once the frame split is known
  std::unique_ptr<FrameCarry> walkCarry;
  auto processFrame = [&](int fr, Cloud &cloud) {
        if (typeFilter <= 0 && cloud.nop == 0) {
          cloud = load(file, cloud.currentFrame, loadType);
        }
        std::ostringstream line;
        FrameCarry fresh(false);
        const int one = runOne(line, walkCarry ? *walkCarry : fresh, cloud);
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
  int lo = frame;
  int hi = last;
  int part = 0;
  int parts = 1;
#ifdef SEAMS_HAS_MPI
  if (ranks.size > 1 && frameSource(file) == FrameSource::dump) {
    part = ranks.rank;
    parts = ranks.size;
  } else if (ranks.rank != 0) {
    // Other readers end a range at its first missing frame, so rank 0
    // reads them alone
    lo = hi + 1;
  }
#endif
  // One worker sees every frame of its range in order; RingUpdater state
  // cannot be shared between workers
  if (carryFrames && parts == 1 &&
      (jobs == 1 || frameSource(file) != FrameSource::dump)) {
    walkCarry = std::make_unique<FrameCarry>(true);
  }
  try {
    forEachInputFrame(file, lo, hi, typeFilter, jobs, processFrame, part, parts);
  } catch (...) {
    std::cerr << colorizer.error(exceptionText()) << "\n";
    return 2;
  }
#ifdef SEAMS_HAS_MPI
  if (ranks.size > 1 &&
      !gatherFrameResults(ranks, outputLines, outputStatus, outputSeen)) {
    if (ranks.rank == 0) {
      std::cerr << colorizer.error("frame results did not gather: past 2 GiB, "
                                   "or a record that does not unpack")
                << "\n";
    }
    return 2;
  }
#endif
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
