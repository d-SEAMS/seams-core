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
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#include <seams_input.hpp>

#include <neighbours.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <iostream>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

#include <sys/stat.h>

#ifdef SEAMS_HAS_CHEMFILES
#include <chemfiles.hpp>
#endif

namespace {

std::string extensionOf(const std::string &path) {
  const auto dot = path.find_last_of('.');
  if (dot == std::string::npos) {
    return {};
  }
  return path.substr(dot + 1);
}

#ifdef SEAMS_HAS_CHEMFILES

bool extensionIs(const std::string &path, const char *ext) {
  const auto got = extensionOf(path);
  const auto n = std::char_traits<char>::length(ext);
  return got.size() == n && std::memcmp(got.data(), ext, n) == 0;
}

bool digitsOnly(const std::string &text) {
  return !text.empty() &&
         text.find_first_not_of("0123456789") == std::string::npos;
}

std::string upperToken(std::string text) {
  for (char &ch : text) {
    ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
  }
  return text;
}

bool oxygenToken(const std::string &token) {
  const std::string name = upperToken(token);
  return name == "O" || name == "OW" || name == "OH";
}

bool hydrogenToken(const std::string &token) {
  const std::string name = upperToken(token);
  if (name == "H" || name == "HW" || name == "HO") {
    return true;
  }
  return name.size() > 2 && name.rfind("HW", 0) == 0 &&
         digitsOnly(name.substr(2));
}

int uffElement(const std::string &token) {
  const auto underscore = token.find('_');
  if (underscore == std::string::npos || underscore < 1 || underscore > 2) {
    return 0;
  }
  std::string symbol = token.substr(0, underscore);
  if (!std::isalpha(static_cast<unsigned char>(symbol[0])) ||
      (symbol.size() == 2 &&
       !std::isalpha(static_cast<unsigned char>(symbol[1])))) {
    return 0;
  }
  symbol[0] =
      static_cast<char>(std::toupper(static_cast<unsigned char>(symbol[0])));
  if (symbol.size() == 2) {
    symbol[1] =
        static_cast<char>(std::tolower(static_cast<unsigned char>(symbol[1])));
  }
  if (const auto number = chemfiles::Atom(symbol).atomic_number()) {
    return static_cast<int>(*number);
  }
  return 0;
}

int mappedSpecies(const chemfiles::Atom &atom,
                  const sinp::ChemfilesLoadOptions &options) {
  for (const auto &token : {atom.name(), atom.type()}) {
    const auto found = options.typeMap.find(token);
    if (found != options.typeMap.end()) {
      return found->second;
    }
  }
  return 0;
}

// Numeric LAMMPS types stay numbers. Element names become atomic numbers,
// and the water names OW/HW do too when the table has no entry for them.
int speciesOf(const chemfiles::Atom &atom,
              const sinp::ChemfilesLoadOptions &options) {
  if (const int mapped = mappedSpecies(atom, options)) {
    return mapped;
  }
  // The token checks run even when the table already has an atomic
  // number. Chemfiles numbers O and H, and also numbers HO as holmium,
  // so a later return would skip the water-name compares entirely.
  const bool oxygenName = oxygenToken(atom.name());
  const bool oxygenType = oxygenToken(atom.type());
  const bool hydrogenName = hydrogenToken(atom.name());
  const bool hydrogenType = hydrogenToken(atom.type());
  // Combine the two names before consulting the atomic number. A PDB
  // atom can carry an element (carbon) and a water atom name (OH), and
  // the number must still win; the combine is what records each side.
  const bool oxygen = oxygenType || oxygenName;
  const bool hydrogen = hydrogenType || hydrogenName;
  if (digitsOnly(atom.type())) {
    return std::stoi(atom.type());
  }
  if (const auto number = atom.atomic_number()) {
    return static_cast<int>(*number);
  }
  if (oxygen) {
    return 8;
  }
  if (hydrogen) {
    return 1;
  }
  if (const int number = uffElement(atom.name())) {
    return number;
  }
  if (const int number = uffElement(atom.type())) {
    return number;
  }
  return 0;
}

int atomicNumberOf(const chemfiles::Atom &atom) {
  if (const auto number = atom.atomic_number()) {
    return static_cast<int>(*number);
  }
  return 0;
}

void rememberId(molSys::PointCloud<molSys::Point<double>, double> &cloud) {
  const int idx = static_cast<int>(cloud.pts.size()) - 1;
  cloud.idIndexMap[cloud.pts.back().atomID] = idx;
}

bool hasName(const std::vector<std::string> &names,
             const std::string &name) {
  return std::find(names.begin(), names.end(), name) != names.end();
}

bool atomWanted(const chemfiles::Topology &topology, std::size_t index,
                const sinp::ChemfilesLoadOptions &options) {
  const std::string &atomName = topology[index].name();
  if ((!options.atomNames.empty() &&
       !hasName(options.atomNames, atomName)) ||
      hasName(options.excludedAtomNames, atomName)) {
    return false;
  }
  std::string residueName;
  if (const auto residue = topology.residue_for_atom(index)) {
    residueName = residue->name();
  }
  if ((!options.residueNames.empty() &&
       !hasName(options.residueNames, residueName)) ||
      hasName(options.excludedResidueNames, residueName)) {
    return false;
  }
  return true;
}

void applyCell(molSys::PointCloud<molSys::Point<double>, double> &cloud,
               const chemfiles::UnitCell &cell) {
  if (cell.shape() == chemfiles::UnitCell::INFINITE) {
    cloud.box.clear();
    cloud.boxLow.clear();
    return;
  }
  const auto lengths = cell.lengths();
  const auto angles = cell.angles();
  nneigh::setBoxFromLengthsAngles(cloud.box, cloud.boxLow, lengths[0],
                                 lengths[1], lengths[2], angles[0], angles[1],
                                 angles[2]);
}

std::string topologyFor(const std::string &filename) {
  const auto ext = extensionOf(filename);
  const auto dot = filename.find_last_of('.');
  std::string topology;
  if (ext == "dcd") {
    topology = filename.substr(0, dot) + ".pdb";
  } else if (ext == "xtc" || ext == "trr" || ext == "nc") {
    topology = filename.substr(0, dot) + ".gro";
  }
  if (!topology.empty()) {
    std::ifstream probe(topology);
    if (probe) {
      return topology;
    }
  }
  return {};
}

#endif

} // namespace

#ifdef SEAMS_HAS_CHEMFILES

// Set when chemfiles rejects the file, so openFrame can try the native reader.
// A selection that keeps no atoms does not set it. The message is the
// chemfiles reason; an XYZ failure is reported from it and is not replaced
// by a native parse of whatever lines happened to look like atoms.
thread_local bool chemfilesUnreadable = false;
thread_local bool framePastEnd = false;
thread_local std::string chemfilesError;

bool extensionFold(const std::string &ext, const char *want) {
  const auto n = std::char_traits<char>::length(want);
  if (ext.size() != n) {
    return false;
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (std::tolower(static_cast<unsigned char>(ext[i])) !=
        std::tolower(static_cast<unsigned char>(want[i]))) {
      return false;
    }
  }
  return true;
}

std::string forcedFormat(const std::string &filename) {
  const auto ext = extensionOf(filename);
  // .cif is small-molecule CIF. .mmcif is mmCIF. Guessing the other
  // way round is how a written mmCIF was opened as CIF.
  if (ext == "mmcif") {
    return "mmCIF";
  }
  if (ext == "cif") {
    return "CIF";
  }
  // LAMMPS Data registers no extension (chemfiles FormatMetadata::extension
  // is empty). The frame count and the frame read both pass this name.
  if (extensionFold(ext, "lammpsdata")) {
    return "LAMMPS Data";
  }
  return "";
}

// chemfiles 0.10 mmCIFFormat::init_ dereferences a null block when the
// file has no data_ record. Reject those files before the constructor.
bool mmcifHasData(const std::string &filename) {
  if (extensionOf(filename) != "mmcif") {
    return true;
  }
  std::ifstream in(filename);
  if (!in) {
    return false;
  }
  std::string chunk(8192, '\0');
  in.read(chunk.data(), static_cast<std::streamsize>(chunk.size()));
  chunk.resize(static_cast<std::size_t>(std::max<std::streamsize>(in.gcount(), 0)));
  return chunk.find("data_") != std::string::npos ||
         chunk.find("_atom_site") != std::string::npos ||
         chunk.find("_cell.") != std::string::npos;
}

#ifdef SEAMS_CHEMFILES_PATCHED
// The vendored chemfiles that Meson builds exports these. The trajectory,
// the frame, and a selection list are created and destroyed there, so
// their free uses the allocator that allocated them.
namespace chemfiles {
int seams_with_frame(const std::string &path, const std::string &format,
                     std::size_t step, const std::string &topology,
                     const std::string &selection, bool read, bool empty_xyz,
                     int (*visit)(const Frame &frame, const std::size_t *ids,
                                  std::size_t nids, void *user),
                     void *user, std::string &error, std::size_t &nsteps);
int seams_each_frame(const std::string &path, const std::string &format,
                     std::size_t first, std::size_t last, const std::string &topology,
                     int (*visit)(const Frame &frame, std::size_t step, void *user),
                     void *user, std::string &error);
}
namespace framevisit = chemfiles;
#else
// An installed chemfiles has no seams_with_frame, so the same visit is
// written here against the public API. cmake/patches/
// chemfiles-0.10.4-frame-visit.patch is the vendored twin; keep the two
// in step. The return codes are the ones the callers switch on.
namespace framevisit {

const char *const mmcifSeqIdMessage =
    "mmCIF atom_site has no label_entity_id column and an unset label_seq_id";

// Splits a CIF data line into tokens. A quote opens a token only at its
// start and closes at the same quote before white space, so a token with a
// space in it stays whole. Quotes are kept, which leaves a quoted "." apart
// from the unset value ".".
void cifTokens(const std::string &line, std::vector<std::string> &out) {
  std::size_t at = 0;
  const std::size_t n = line.size();
  while (at < n) {
    while (at < n && std::isspace(static_cast<unsigned char>(line[at]))) {
      ++at;
    }
    if (at >= n) {
      break;
    }
    std::size_t end = at + 1;
    if (line[at] == '\'' || line[at] == '"') {
      const char quote = line[at];
      while (end < n &&
             !(line[end] == quote && (end + 1 == n || std::isspace(static_cast<unsigned char>(line[end + 1]))))) {
        ++end;
      }
      end = std::min(end + 1, n);
    } else {
      while (end < n && !std::isspace(static_cast<unsigned char>(line[end]))) {
        ++end;
      }
    }
    out.emplace_back(line.substr(at, end - at));
    at = end;
  }
}

// chemfiles 0.10.4 mmCIFFormat::read dereferences the label_entity_id
// column to classify a residue. When the atom_site loop has no such column
// and an atom has no label_seq_id, which a CIF writes as "." or "?", the
// dereference is of a missing column and the process dies. chemfiles 0.10.4
// writes exactly that for an atom without a residue, so its own output is
// refused here. The vendored copy reads such a file
// (cmake/patches/chemfiles-0.10.4-mmcif.patch).
bool mmcifSeqIdUnset(const std::string &path, const std::string &format) {
  if (upperToken(format) != "MMCIF" &&
      !(format.empty() && extensionIs(path, "mmcif"))) {
    return false;
  }
  std::ifstream in(path);
  std::string line;
  std::vector<std::string> columns;
  std::vector<std::string> tokens;
  bool header = false;
  bool rows = false;
  std::size_t seqColumn = 0;
  while (std::getline(in, line)) {
    const auto first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos) {
      continue;
    }
    std::string text = line.substr(first);
    std::string lowered = text;
    for (char &ch : lowered) {
      ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }
    if (rows) {
      // The loop ends at the next keyword, tag or comment line.
      if (text[0] == '_' || text[0] == '#' || lowered.rfind("loop_", 0) == 0 ||
          lowered.rfind("data_", 0) == 0) {
        return false;
      }
      cifTokens(text, tokens);
      while (tokens.size() >= columns.size()) {
        const std::string &value = tokens[seqColumn];
        if (value == "." || value == "?") {
          return true;
        }
        tokens.erase(tokens.begin(),
                     tokens.begin() + static_cast<std::ptrdiff_t>(columns.size()));
      }
      continue;
    }
    if (lowered.rfind("loop_", 0) == 0) {
      header = true;
      columns.clear();
      continue;
    }
    if (header && text[0] == '_') {
      const auto stop = lowered.find_first_of(" \t\r");
      columns.emplace_back(lowered.substr(0, stop));
      continue;
    }
    if (header) {
      header = false;
      if (columns.empty() || columns.front().rfind("_atom_site.", 0) != 0) {
        continue;
      }
      std::size_t seq = columns.size();
      for (std::size_t i = 0; i < columns.size(); ++i) {
        if (columns[i] == "_atom_site.label_entity_id") {
          return false;
        }
        if (columns[i] == "_atom_site.label_seq_id") {
          seq = i;
        }
      }
      if (seq == columns.size()) {
        return false;
      }
      seqColumn = seq;
      rows = true;
      tokens.clear();
      cifTokens(text, tokens);
      while (tokens.size() >= columns.size()) {
        if (tokens[seqColumn] == "." || tokens[seqColumn] == "?") {
          return true;
        }
        tokens.erase(tokens.begin(),
                     tokens.begin() + static_cast<std::ptrdiff_t>(columns.size()));
      }
    }
  }
  return false;
}

// chemfiles 0.10.4 reads counts[0] of an empty split when a MOL2 molecule
// has no counts line, and segfaults. The vendored copy throws instead
// (cmake/patches/chemfiles-0.10.4-mol2.patch). An installed copy is not
// probed when cross compiling, so the same files are refused here.
bool mol2CountsMissing(const std::string &path, const std::string &format) {
  if (upperToken(format) != "MOL2" &&
      !(format.empty() && extensionIs(path, "mol2"))) {
    return false;
  }
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    if (line.rfind("@<TRIPOS>MOLECULE", 0) != 0) {
      continue;
    }
    std::getline(in, line);
    if (!std::getline(in, line) ||
        line.find_first_not_of(" \r") == std::string::npos) {
      return true;
    }
  }
  return false;
}

int seams_with_frame(const std::string &path, const std::string &format,
                     std::size_t step, const std::string &topology,
                     const std::string &selection, bool read, bool empty_xyz,
                     int (*visit)(const chemfiles::Frame &frame,
                                  const std::size_t *ids, std::size_t nids,
                                  void *user),
                     void *user, std::string &error, std::size_t &nsteps) {
  try {
    if (mol2CountsMissing(path, format)) {
      error = "missing the number of atoms for MOL2";
      return 4;
    }
    if (mmcifSeqIdUnset(path, format)) {
      error = mmcifSeqIdMessage;
      return 4;
    }
    chemfiles::Trajectory trajectory(path, 'r', format);
    nsteps = trajectory.nsteps();
    if (!read) {
      return 0;
    }
    if (empty_xyz && nsteps == 0) {
      error = "XYZ read failed: the file has no frames";
      return 6;
    }
    if (nsteps > 0 && step >= nsteps) {
      return 2;
    }
    if (!topology.empty()) {
      trajectory.set_topology(topology);
    }
    chemfiles::Frame frame = trajectory.read_step(step);
    const std::size_t *ids = nullptr;
    std::size_t nids = static_cast<std::size_t>(-1);
    std::vector<std::size_t> picked;
    if (!selection.empty()) {
      chemfiles::Selection chosen(selection);
      picked = chosen.list(frame);
      ids = picked.data();
      nids = picked.size();
    }
    if (visit == nullptr) {
      return 0;
    }
    return visit(frame, ids, nids, user);
  } catch (const chemfiles::SelectionError &err) {
    error = err.what();
    return 3;
  } catch (const chemfiles::Error &err) {
    error = err.what();
    return 4;
  } catch (const std::exception &err) {
    error = err.what();
    return 5;
  }
}

int seams_each_frame(const std::string &path, const std::string &format,
                     std::size_t first, std::size_t last,
                     const std::string &topology,
                     int (*visit)(const chemfiles::Frame &frame,
                                  std::size_t step, void *user),
                     void *user, std::string &error) {
  try {
    if (mol2CountsMissing(path, format)) {
      error = "missing the number of atoms for MOL2";
      return 4;
    }
    if (mmcifSeqIdUnset(path, format)) {
      error = mmcifSeqIdMessage;
      return 4;
    }
    chemfiles::Trajectory trajectory(path, 'r', format);
    if (!topology.empty()) {
      trajectory.set_topology(topology);
    }
    const std::size_t nsteps = trajectory.nsteps();
    for (std::size_t step = first; step < nsteps && step <= last; ++step) {
      // read_step for every step: a format's own read position need not follow
      // a read_step, so read() after it can return another step.
      chemfiles::Frame frame = trajectory.read_step(step);
      if (visit != nullptr && visit(frame, step, user) != 0) {
        break;
      }
    }
    return 0;
  } catch (const chemfiles::Error &err) {
    error = err.what();
    return 4;
  } catch (const std::exception &err) {
    error = err.what();
    return 5;
  }
}

} // namespace framevisit
#endif

int sinp::chemfilesSteps(const std::string &filename) {
  if (!mmcifHasData(filename)) {
    chemfilesError = "mmCIF file has no data_ block";
    return -1;
  }
  std::string error;
  std::size_t nsteps = 0;
  const int rc = framevisit::seams_with_frame(
      filename, forcedFormat(filename), 0, "", "", false, false, nullptr,
      nullptr, error, nsteps);
  if (rc != 0) {
    chemfilesError = error.empty() ? currentExceptionText() : error;
    return -1;
  }
  return static_cast<int>(nsteps);
}

molSys::PointCloud<molSys::Point<double>, double>
sinp::readChemfiles(std::string filename, int targetFrame,
                    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                    int typeFilter) {
  ChemfilesRead request;
  request.frame = targetFrame;
  request.typeFilter = typeFilter;
  return readChemfiles(std::move(filename), yCloud, request);
}

struct ChemVisit {
  const sinp::ChemfilesRead *request;
  molSys::PointCloud<molSys::Point<double>, double> *cloud;
};

int visitChemFrame(const chemfiles::Frame &frame, const std::size_t *ids,
                   std::size_t nids, void *user) {
  auto *ctx = static_cast<ChemVisit *>(user);
  const sinp::ChemfilesRead &request = *ctx->request;
  auto &yCloud = *ctx->cloud;
  applyCell(yCloud, frame.cell());
    const auto &topology = frame.topology();
    const auto positions = frame.positions();
    const std::size_t n = frame.size();

    std::vector<char> wanted(n, 1);
    if (nids != static_cast<std::size_t>(-1)) {
      std::fill(wanted.begin(), wanted.end(), 0);
      for (std::size_t k = 0; k < nids; ++k) {
        if (ids[k] < n) {
          wanted[ids[k]] = 1;
        }
      }
    }

    bool anyOxygen = false;
    bool anyElement = false;
    int numericTwos = 0;
    for (std::size_t i = 0; i < n; ++i) {
      const int number = atomicNumberOf(topology[i]);
      if (!atomWanted(topology, i, request.load)) {
        wanted[i] = 0;
        continue;
      }
      const int species = speciesOf(topology[i], request.load);
      const bool oxygenName = oxygenToken(topology[i].name());
      const bool oxygenType = oxygenToken(topology[i].type());
      const bool hydrogenName = hydrogenToken(topology[i].name());
      const bool hydrogenType = hydrogenToken(topology[i].type());
      const bool named = ((number > 0) | oxygenName | oxygenType | hydrogenName |
                          hydrogenType) != 0;
      if (named) {
        anyElement = true;
      }
      if (number == 8 || species == 8) {
        anyOxygen = true;
      }
      if (!named && species == 2) {
        numericTwos += 1;
      }
    }
    const bool keepOxygen = request.guessOxygen && anyOxygen;
    const bool keepTwo =
        request.guessOxygen && !anyElement && numericTwos > 0;
    const bool keepOne =
        request.guessOxygen && !anyElement && numericTwos == 0;

    for (std::size_t i = 0; i < n; ++i) {
      if (!wanted[i]) {
        continue;
      }
      int species = speciesOf(topology[i], request.load);
      const int number = atomicNumberOf(topology[i]);
      // A file with no names and no elements is one species.
      // A GROMACS TRR is that case, and the live module labels it type 1.
      if (!anyElement && species == 0) {
        species = 1;
      }
      if (keepOxygen && number != 8 && species != 8) {
        continue;
      }
      if (keepTwo && species != 2) {
        continue;
      }
      if (keepOne && species != 1) {
        continue;
      }
      if (request.typeFilter > 0 && species != request.typeFilter) {
        continue;
      }
      molSys::Point<double> point;
      point.type = species;
      point.atomID = static_cast<int>(i) + 1;
      point.molID = point.atomID;
      // LAMMPS writes mol as a custom column. Chemfiles keeps it on the atom.
      if (auto mol = topology[i].get<chemfiles::Property::DOUBLE>("mol")) {
        point.molID = static_cast<int>(std::lround(*mol));
      } else if (const auto residue = topology.residue_for_atom(i)) {
        if (const auto residueId = residue->id()) {
          point.molID = static_cast<int>(*residueId);
        }
      }
      point.x = positions[i][0];
      point.y = positions[i][1];
      point.z = positions[i][2];
      point.inSlice = true;
      yCloud.pts.push_back(point);
      rememberId(yCloud);
    }
    yCloud.nop = static_cast<int>(yCloud.pts.size());
    yCloud.currentFrame = request.frame;
    return 0;
}

molSys::PointCloud<molSys::Point<double>, double>
sinp::readChemfiles(std::string filename,
                    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                    const ChemfilesRead &request) {
  chemfilesUnreadable = false;
  chemfilesError.clear();
  yCloud = molSys::clearPointCloud(yCloud);
  if (request.frame < 1) {
    return yCloud;
  }
  if (!mmcifHasData(filename)) {
    chemfilesUnreadable = true;
    chemfilesError = "mmCIF file has no data_ block";
    return yCloud;
  }
  const std::string topology = topologyFor(filename);
  ChemVisit ctx{&request, &yCloud};
  std::string error;
  std::size_t nsteps = 0;
  const int rc = framevisit::seams_with_frame(
      filename, forcedFormat(filename),
      static_cast<std::size_t>(request.frame - 1), topology, request.selection,
      true, extensionIs(filename, "xyz"), visitChemFrame, &ctx, error, nsteps);
  if (rc == 6) {
    throw std::runtime_error(error);
  }
  if (rc == 2) {
    framePastEnd = true;
    yCloud = molSys::clearPointCloud(yCloud);
    return yCloud;
  }
  if (rc == 3) {
    std::cerr << "chemfiles selection on " << filename << ": " << error << "\n";
    yCloud = molSys::clearPointCloud(yCloud);
    return yCloud;
  }
  if (rc != 0) {
    chemfilesUnreadable = true;
    chemfilesError = error;
    yCloud = molSys::clearPointCloud(yCloud);
    return yCloud;
  }
  return yCloud;
}

struct ChemRange {
  sinp::ChemfilesRead request;
  const std::function<void(int, molSys::PointCloud<molSys::Point<double>, double> &)> *fn;
  int visited = 0;
};

int visitChemRange(const chemfiles::Frame &frame, std::size_t step, void *user) {
  auto *range = static_cast<ChemRange *>(user);
  range->request.frame = static_cast<int>(step) + 1;
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  ChemVisit ctx{&range->request, &cloud};
  visitChemFrame(frame, nullptr, static_cast<std::size_t>(-1), &ctx);
  range->visited += 1;
  (*range->fn)(range->request.frame, cloud);
  return 0;
}

void sinp::forEachChemfilesStep(
    const std::string &filename, int first, int last, int typeFilter,
    const ChemfilesLoadOptions &options,
    const std::function<void(
        int, molSys::PointCloud<molSys::Point<double>, double> &)> &fn) {
  const std::string topology = topologyFor(filename);
  ChemRange range;
  range.request.typeFilter = typeFilter;
  range.request.load = options;
  range.fn = &fn;
  const std::size_t from = static_cast<std::size_t>(std::max(first, 1) - 1);
  const std::size_t to = last > 0 ? static_cast<std::size_t>(last - 1)
                                  : std::numeric_limits<std::size_t>::max();
  std::string error;
  int rc = 0;
  if (mmcifHasData(filename)) {
    rc = framevisit::seams_each_frame(filename, forcedFormat(filename), from,
                                      to, topology, visitChemRange, &range,
                                      error);
  } else {
    rc = 4;
  }
  // An unreadable file gives the first frame as an empty cloud, as one
  // readChemfiles call per frame did.
  if (rc != 0 && range.visited == 0) {
    molSys::PointCloud<molSys::Point<double>, double> empty;
    fn(std::max(first, 1), empty);
  }
}

void sinp::forEachChemfilesStep(
    const std::string &filename, int first, int last, int typeFilter,
    const std::function<void(
        int, molSys::PointCloud<molSys::Point<double>, double> &)> &fn) {
  forEachChemfilesStep(filename, first, last, typeFilter,
                       ChemfilesLoadOptions{}, fn);
}

bool sinp::forEachChemfilesFrame(
    const std::string &filename, int first, int last, int typeI,
    const std::string &selection, int part, int parts,
    const std::function<void(
        int, molSys::PointCloud<molSys::Point<double>, double> &)> &fn) {
  if (extensionOf(filename) == "con") {
    return false;
  }
  const int steps = chemfilesSteps(filename);
  if (steps < 0) {
    return false;
  }
  if (steps == 0 && extensionIs(filename, "xyz")) {
    throw std::runtime_error("XYZ read failed: the file has no frames");
  }
  const int begin = std::max(1, first);
  int end = last;
  if (end <= 0 || end > steps) {
    end = steps;
  }
  if (steps > 0 && last > 0 && begin > end) {
    throw std::runtime_error(
        "frame " + std::to_string(first) +
        " is past the last frame in this file");
  }
  if (parts < 1) {
    parts = 1;
  }
  if (part < 0 || part >= parts) {
    return true;
  }
  for (int frame = begin; frame <= end; ++frame) {
    if (((frame - begin) % parts) != part) {
      continue;
    }
    auto cloud = openFrame(filename, frame, typeI, selection);
    fn(frame, cloud);
  }
  return true;
}

#endif

std::string sinp::currentExceptionText() {
  try {
    throw;
  } catch (const std::exception &err) {
    return err.what();
  } catch (...) {
    return "unhandled error";
  }
}

#ifdef SEAMS_HAS_CHEMFILES
namespace {

// One bit per frame: ids are exactly {1..N} and the box has no tilt.
// Chemfiles then stores id k at index k-1. A gap would index past the
// atom count. Chemfiles 0.10.4 also takes xhi_bound - xlo_bound as the
// edge of a tilted cell without removing the tilt padding, so those
// frames stay on the native reader too. The file is walked once and kept.

bool startsWith(const char *line, const char *prefix) {
  const auto n = std::char_traits<char>::length(prefix);
  return std::strncmp(line, prefix, n) == 0;
}

struct DumpIdIndex {
  off_t bytes = -1;
  std::filesystem::file_time_type mtime{};
  std::vector<char> oneToN;
};

std::mutex gDumpIdMu;
std::unordered_map<std::string, DumpIdIndex> gDumpIds;

DumpIdIndex scanDumpIds(const std::string &path) {
  DumpIdIndex out;
  struct stat st {};
  if (stat(path.c_str(), &st) != 0) {
    return out;
  }
  out.bytes = st.st_size;
  std::error_code ec;
  out.mtime = std::filesystem::last_write_time(path, ec);
  // getline and operator>> allocate inside libstdc++. A host that loads
  // this library after that runtime (Julia) then frees those buffers
  // with the sanitizer allocator. Read and split the line here instead.
  std::FILE *in = std::fopen(path.c_str(), "r");
  if (in == nullptr) {
    return out;
  }
  const auto t0 = std::chrono::steady_clock::now();
  int seen = 0;
  int natoms = 0;
  int idCol = -1;
  int got = 0;
  bool atoms = false;
  bool tilted = false;
  // Ids in [1, natoms], one per row read. The header count is untrusted,
  // so the bitmap is sized only once the rows have shown it is real.
  std::vector<int> ids;
  std::vector<char> seenId;
  // The characters live in a buffer this function mallocs and frees.
  // A std::string grown by iostream would be allocated in libstdc++ and
  // freed by the sanitizer build of this function.
  struct DumpLine {
    char *data = nullptr;
    std::size_t len = 0;
    std::size_t cap = 0;
    void clear() {
      len = 0;
      if (data != nullptr) {
        data[0] = '\0';
      }
    }
    void push(char ch) {
      if (len + 1 >= cap) {
        const std::size_t next = cap == 0 ? 128 : cap * 2;
        char *grown = static_cast<char *>(std::realloc(data, next));
        if (grown == nullptr) {
          return;
        }
        data = grown;
        cap = next;
      }
      data[len++] = ch;
      data[len] = '\0';
    }
    const char *c_str() const { return data != nullptr ? data : ""; }
    ~DumpLine() { std::free(data); }
    DumpLine() = default;
    DumpLine(const DumpLine &) = delete;
    DumpLine &operator=(const DumpLine &) = delete;
  } line;
  auto finish = [&]() {
    if (seen <= 0) {
      return;
    }
    out.oneToN.resize(static_cast<std::size_t>(seen) + 1, 0);
    bool ok = natoms > 0 && ids.size() == static_cast<std::size_t>(natoms) && !tilted;
    if (ok) {
      seenId.assign(static_cast<std::size_t>(natoms) + 1, 0);
      for (const int id : ids) {
        if (seenId[static_cast<std::size_t>(id)] != 0) {
          ok = false;
          break;
        }
        seenId[static_cast<std::size_t>(id)] = 1;
      }
    }
    out.oneToN[static_cast<std::size_t>(seen)] = ok ? 1 : 0;
  };
  auto nextLine = [&]() {
    line.clear();
    int ch = std::fgetc(in);
    if (ch == EOF) {
      return false;
    }
    while (ch != EOF && ch != '\n') {
      line.push(static_cast<char>(ch));
      ch = std::fgetc(in);
    }
    return true;
  };
  while (nextLine()) {
    if (startsWith(line.c_str(), "ITEM: TIMESTEP")) {
      finish();
      seen += 1;
      atoms = false;
      natoms = 0;
      idCol = -1;
      got = 0;
      ids.clear();
      tilted = false;
      seenId.clear();
      continue;
    }
    if (startsWith(line.c_str(), "ITEM: BOX BOUNDS")) {
      if (std::strstr(line.c_str(), "xy") != nullptr) {
        for (int axis = 0; axis < 3 && nextLine(); ++axis) {
          double lo = 0.0;
          double hi = 0.0;
          double tilt = 0.0;
          if (std::sscanf(line.c_str(), "%lf %lf %lf", &lo, &hi, &tilt) == 3 && tilt != 0.0) {
            tilted = true;
          }
        }
      }
      continue;
    }
    if (startsWith(line.c_str(), "ITEM: NUMBER OF ATOMS")) {
      if (nextLine()) {
        natoms = std::atoi(line.c_str());
      }
      continue;
    }
    if (startsWith(line.c_str(), "ITEM: ATOMS")) {
      idCol = -1;
      const char *cursor = line.c_str();
      auto nextField = [&](const char *&p, const char *&start) {
        while (*p == ' ' || *p == '\t' || *p == '\r') {
          ++p;
        }
        start = p;
        while (*p != '\0' && *p != ' ' && *p != '\t' && *p != '\r') {
          ++p;
        }
        return p != start;
      };
      const char *field = nullptr;
      (void)nextField(cursor, field);
      (void)nextField(cursor, field);
      int col = 0;
      while (nextField(cursor, field)) {
        if (cursor - field == 2 && field[0] == 'i' && field[1] == 'd') {
          idCol = col;
        }
        col += 1;
      }
      atoms = true;
      got = 0;
      ids.clear();
      continue;
    }
    if (!atoms || got >= natoms) {
      continue;
    }
    const char *field = line.c_str();
    int fieldCount = 0;
    const char *idStart = nullptr;
    const char *idEnd = nullptr;
    while (*field != '\0') {
      while (*field == ' ' || *field == '\t' || *field == '\r') {
        ++field;
      }
      if (*field == '\0') {
        break;
      }
      const char *start = field;
      while (*field != '\0' && *field != ' ' && *field != '\t' && *field != '\r') {
        ++field;
      }
      if (fieldCount == idCol) {
        idStart = start;
        idEnd = field;
      }
      fieldCount += 1;
    }
    got += 1;
    if (idStart == nullptr) {
      continue;
    }
    char idText[32];
    const std::size_t idChars = std::min(sizeof(idText) - 1, static_cast<std::size_t>(idEnd - idStart));
    std::memcpy(idText, idStart, idChars);
    idText[idChars] = '\0';
    const int id = std::atoi(idText);
    if (id >= 1 && id <= natoms) {
      ids.push_back(id);
    }
  }
  finish();
  std::fclose(in);
  const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                      std::chrono::steady_clock::now() - t0)
                      .count();
  if (const char *log = std::getenv("SEAMS_MOL_INDEX_LOG")) {
    if (log[0] != '\0' && log[0] != '0') {
      std::cerr << "mol index frames " << seen << " scan_ms " << ms << "\n";
    }
  }
  return out;
}

bool dumpFrameIsOneToN(const std::string &path, int frame) {
  if (frame < 1) {
    return false;
  }
  struct stat st {};
  const bool exists = stat(path.c_str(), &st) == 0;
  std::error_code ec;
  const auto mtime = std::filesystem::last_write_time(path, ec);
  std::lock_guard<std::mutex> lock(gDumpIdMu);
  auto it = gDumpIds.find(path);
  const bool fresh = exists && it != gDumpIds.end() && it->second.bytes == st.st_size &&
                     !ec && it->second.mtime == mtime;
  if (!fresh) {
    gDumpIds[path] = scanDumpIds(path);
    it = gDumpIds.find(path);
  }
  if (frame >= static_cast<int>(it->second.oneToN.size())) {
    return false;
  }
  return it->second.oneToN[static_cast<std::size_t>(frame)] != 0;
}

} // namespace
#endif

molSys::PointCloud<molSys::Point<double>, double>
sinp::openFrame(const std::string &filename, int frame, int typeI,
                const std::string &selection,
                const ChemfilesLoadOptions &options) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  if (extensionOf(filename) == "con") {
#ifdef SEAMS_HAS_READCON
    return readCon(filename, frame, cloud);
#else
    throw std::runtime_error(
        "readcon-core is not in this build; CON input needs it");
#endif
  }
#ifdef SEAMS_HAS_CHEMFILES
  ChemfilesRead request;
  request.frame = frame;
  request.selection = selection;
  request.load = options;
  if (!selection.empty()) {
    request.typeFilter = typeI > 0 ? typeI : -1;
  } else if (typeI < 0) {
    request.typeFilter = -1;
  } else if (typeI > 0) {
    request.typeFilter = typeI;
  } else {
    request.guessOxygen = true;
  }
  const auto extNow = extensionOf(filename);
  const bool dumpFile = extNow == "lammpstrj" || extNow == "dump";
  if (dumpFile && !dumpFrameIsOneToN(filename, frame)) {
      // Gaps, ids that do not start at 1, or a tilted box (see
      // scanDumpIds). One native read.
      if (typeI < 0) {
        return readLammpsTrj(filename, frame, cloud);
      }
      if (typeI > 0) {
        return readLammpsTrjO(filename, frame, cloud, typeI);
      }
      cloud = readLammpsTrjO(filename, frame, cloud, 2);
      if (cloud.nop == 0) {
        molSys::PointCloud<molSys::Point<double>, double> again;
        cloud = readLammpsTrjO(filename, frame, again, 1);
      }
      return cloud;
  }
  framePastEnd = false;
  cloud = readChemfiles(filename, cloud, request);
  if (framePastEnd) {
    throw std::runtime_error(
        "frame " + std::to_string(frame) +
        " is past the last frame in this file");
  }
  if (!chemfilesUnreadable) {
    if (dumpFile) {
      // Chemfiles reports the bound span as lattice length a and the
      // tilt as the cell angle. setBoxFromLengthsAngles then treats
      // that span as an edge and pads it again, so the recovered lx
      // is the span (31.05 A on genice sH, not the 24.84 A edge).
      // The native BOX BOUNDS parser stores the span and the tilt once.
      try {
        molSys::PointCloud<molSys::Point<double>, double> native;
        native = readLammpsTrj(filename, frame, native);
        if (native.box.size() >= 3) {
          cloud.box = std::move(native.box);
          cloud.boxLow = std::move(native.boxLow);
        }
      } catch (const std::exception &) {
      }
    }
    return cloud;
  }
  {
    const auto ext = extensionOf(filename);
    if (ext == "xyz") {
      // Chemfiles rejected the file. The native reader is one frame.
      if (frame > 1) {
        throw std::runtime_error(
            "frame " + std::to_string(frame) +
            " is past the last frame in this file");
      }
      return readXYZ(filename);
    }
    if (ext == "lammpstrj" || ext == "dump") {
      if (typeI < 0) {
        return readLammpsTrj(filename, frame, cloud);
      }
      if (typeI > 0) {
        return readLammpsTrjO(filename, frame, cloud, typeI);
      }
      cloud = readLammpsTrjO(filename, frame, cloud, 2);
      if (cloud.nop == 0) {
        molSys::PointCloud<molSys::Point<double>, double> again;
        cloud = readLammpsTrjO(filename, frame, again, 1);
      }
      return cloud;
    }
    if (chemfilesError.empty()) {
      throw std::runtime_error("chemfiles cannot read " + filename);
    }
    throw std::runtime_error(chemfilesError);
  }
#else
  (void)selection;
  const auto ext = extensionOf(filename);
  std::string lowerExt = ext;
  for (char &ch : lowerExt) {
    ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
  }
  // These formats are read only through chemfiles. Without it they are a
  // refused read, as a chemfiles build reports them, not a LAMMPS dump
  // with no frames.
  static const char *const chemfilesOnly[] = {
      "mol2", "mmcif", "cif", "pdb", "gro", "dcd", "xtc", "trr", "tng", "nc",
      "netcdf", "sdf", "cml", "cssr", "arc", "mmtf", "lammpsdata"};
  for (const char *only : chemfilesOnly) {
    if (lowerExt == only) {
      throw std::runtime_error("chemfiles is not in this build; " + lowerExt +
                               " input needs it");
    }
  }
  if (ext == "xyz") {
    if (frame > 1) {
      throw std::runtime_error(
          "frame " + std::to_string(frame) +
          " is past the last frame in this file");
    }
    return readXYZ(filename);
  }
  if (typeI < 0) {
    return readLammpsTrj(filename, frame, cloud);
  }
  if (typeI > 0) {
    return readLammpsTrjO(filename, frame, cloud, typeI);
  }
  cloud = readLammpsTrjO(filename, frame, cloud, 2);
  if (cloud.nop > 0) {
    return cloud;
  }
  molSys::PointCloud<molSys::Point<double>, double> again;
  return readLammpsTrjO(filename, frame, again, 1);
#endif
}

#ifndef SEAMS_HAS_CHEMFILES
int sinp::chemfilesSteps(const std::string &) { return -1; }
#endif
