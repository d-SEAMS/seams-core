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
#include <cstdlib>
#include <cstring>
#include <fstream>
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

bool extensionIs(const std::string &path, const char *ext) {
  const auto got = extensionOf(path);
  const auto n = std::char_traits<char>::length(ext);
  return got.size() == n && std::memcmp(got.data(), ext, n) == 0;
}

#ifdef SEAMS_HAS_CHEMFILES

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
  return name == "H" || name == "HW" || name == "HO";
}

// Numeric LAMMPS types stay numbers. Element names become atomic numbers,
// and the water names OW/HW do too when the table has no entry for them.
int speciesOf(const chemfiles::Atom &atom) {
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
  return "";
}

int sinp::chemfilesSteps(const std::string &filename) {
  try {
    chemfiles::Trajectory trajectory(filename, 'r', forcedFormat(filename));
    return static_cast<int>(trajectory.nsteps());
  } catch (...) {
    chemfilesError = currentExceptionText();
    return -1;
  }
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
  try {
    chemfiles::Trajectory trajectory(filename, 'r', forcedFormat(filename));
    // A DCD has coordinates and a cell, not atom names. A sibling PDB
    // carries the names the host used to choose the analysed species.
    if (extensionIs(filename, "xyz") && trajectory.nsteps() == 0) {
      throw std::runtime_error("XYZ read failed: the file has no frames");
    }
    if (extensionOf(filename) == "dcd") {
      // The suffix check already found a dot.
      const auto dot = filename.find_last_of('.');
      const auto pdb = filename.substr(0, dot) + ".pdb";
      std::ifstream probe(pdb);
      if (probe) {
        trajectory.set_topology(pdb);
      }
    }
    if (trajectory.nsteps() > 0 &&
        static_cast<std::size_t>(request.frame) > trajectory.nsteps()) {
      framePastEnd = true;
      yCloud = molSys::clearPointCloud(yCloud);
      return yCloud;
    }
    auto frame =
        trajectory.read_step(static_cast<std::size_t>(request.frame - 1));
    applyCell(yCloud, frame.cell());
    const auto &topology = frame.topology();
    const auto positions = frame.positions();
    const std::size_t n = frame.size();

    std::vector<char> wanted(n, 1);
    if (!request.selection.empty()) {
      std::fill(wanted.begin(), wanted.end(), 0);
      const chemfiles::Selection selection(request.selection);
      for (const std::size_t index : selection.list(frame)) {
        if (index < n) {
          wanted[index] = 1;
        }
      }
    }

    bool anyOxygen = false;
    bool anyElement = false;
    int numericTwos = 0;
    for (std::size_t i = 0; i < n; ++i) {
      const int number = atomicNumberOf(topology[i]);
      const int species = speciesOf(topology[i]);
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
      int species = speciesOf(topology[i]);
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
    return yCloud;
  } catch (const chemfiles::SelectionError &error) {
    std::cerr << "chemfiles selection on " << filename << ": " << error.what()
              << "\n";
    yCloud = molSys::clearPointCloud(yCloud);
    return yCloud;
  } catch (const chemfiles::Error &error) {
    chemfilesUnreadable = true;
    chemfilesError = error.what();
    yCloud = molSys::clearPointCloud(yCloud);
    return yCloud;
  }
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

namespace {

// One bit per frame: ids are exactly {1..N}. Chemfiles then stores
// id k at index k-1. A gap would index past the atom count, so those
// frames stay on the native reader. The file is walked once and kept.

bool startsWith(const std::string &line, const char *prefix) {
  const auto n = std::char_traits<char>::length(prefix);
  return line.size() >= n && line.compare(0, n, prefix) == 0;
}

struct DumpIdIndex {
  off_t bytes = -1;
  timespec mtime {};
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
  out.mtime = st.st_mtim;
  std::ifstream in(path);
  const auto t0 = std::chrono::steady_clock::now();
  int seen = 0;
  int natoms = 0;
  int idCol = -1;
  int got = 0;
  int unique = 0;
  bool atoms = false;
  std::vector<char> seenId;
  std::string line;
  auto finish = [&]() {
    if (seen <= 0) {
      return;
    }
    out.oneToN.resize(static_cast<std::size_t>(seen) + 1, 0);
    const bool ok = natoms > 0 && unique == natoms;
    out.oneToN[static_cast<std::size_t>(seen)] = ok ? 1 : 0;
  };
  while (std::getline(in, line)) {
    if (startsWith(line, "ITEM: TIMESTEP")) {
      finish();
      seen += 1;
      atoms = false;
      natoms = 0;
      idCol = -1;
      got = 0;
      unique = 0;
      seenId.clear();
      continue;
    }
    if (startsWith(line, "ITEM: NUMBER OF ATOMS")) {
      if (std::getline(in, line)) {
        natoms = std::atoi(line.c_str());
      }
      continue;
    }
    if (startsWith(line, "ITEM: ATOMS")) {
      idCol = -1;
      std::istringstream head(line);
      std::string tok;
      head >> tok >> tok;
      int col = 0;
      while (head >> tok) {
        if (tok == "id") {
          idCol = col;
        }
        col += 1;
      }
      atoms = true;
      got = 0;
      unique = 0;
      seenId.assign(static_cast<std::size_t>(std::max(natoms, 0)) + 1, 0);
      continue;
    }
    if (!atoms || got >= natoms) {
      continue;
    }
    std::istringstream row(line);
    std::vector<std::string> fields;
    std::string tok;
    while (row >> tok) {
      fields.push_back(tok);
    }
    got += 1;
    if (idCol < 0 || idCol >= static_cast<int>(fields.size())) {
      continue;
    }
    const int id = std::atoi(fields[static_cast<std::size_t>(idCol)].c_str());
    if (id >= 1 && id <= natoms && seenId[static_cast<std::size_t>(id)] == 0) {
      seenId[static_cast<std::size_t>(id)] = 1;
      unique += 1;
    }
  }
  finish();
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
  std::lock_guard<std::mutex> lock(gDumpIdMu);
  auto it = gDumpIds.find(path);
  const bool fresh = exists && it != gDumpIds.end() && it->second.bytes == st.st_size &&
                     it->second.mtime.tv_sec == st.st_mtim.tv_sec &&
                     it->second.mtime.tv_nsec == st.st_mtim.tv_nsec;
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

molSys::PointCloud<molSys::Point<double>, double>
sinp::openFrame(const std::string &filename, int frame, int typeI,
                const std::string &selection) {
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
      // Gaps or ids that do not start at 1. Chemfiles stores an atom at
      // index id-1 and would address past the atom count. One native read.
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
    std::cerr << "chemfiles cannot read " << filename << "\n";
    return molSys::clearPointCloud(cloud);
  }
#else
  (void)selection;
  const auto ext = extensionOf(filename);
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
