//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <generic.hpp>
#include <memory>
#include <mutex>
#include <seams_input.hpp>
#include <unordered_map>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {
/**
 * @details Record the ID-to-index entry for the most recently appended point.
 *  A duplicate atom ID means the input is corrupt (two atoms cannot share an
 *  ID within one frame); the map can only keep one of the colliding atoms, so
 *  name the ID loudly instead of failing silently downstream.
 */
void mapAtomIdToIndex(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  const int idx = static_cast<int>(yCloud.pts.size()) - 1;
  const int id = yCloud.pts.back().atomID;
  if (!yCloud.idIndexMap.emplace(id, idx).second) {
    std::cerr << "Warning: duplicate atom ID " << id
              << " in the input frame; the ID-to-index map keeps only the "
                 "last atom read with this ID\n";
    yCloud.idIndexMap[id] = idx;
  }
}

bool isLammpsTimestep(const std::string &line) {
  return line.size() >= 14 && line.compare(0, 14, "ITEM: TIMESTEP") == 0;
}

// Live dump cursor, matching LAMMPS ReaderNative (rerun keeps FILE* at
// the next snapshot) and chemfiles LAMMPSTrajectory::read_next. Random
// access seeks a cached ITEM: TIMESTEP offset table.
struct LammpsDumpSession {
  std::mutex mu;
  std::string path;
  std::ifstream file;
  std::vector<std::uint64_t> offsets;
  int lastFrame = 0;
  std::filesystem::file_time_type mtime{};
  std::uintmax_t size{0};

  bool open(const std::string &filename) {
    path = filename;
    std::error_code ec;
    size = std::filesystem::file_size(filename, ec);
    if (ec) {
      return false;
    }
    mtime = std::filesystem::last_write_time(filename, ec);
    if (ec) {
      return false;
    }
    file.open(filename, std::ios::in | std::ios::binary);
    offsets.clear();
    lastFrame = 0;
    return file.is_open();
  }

  bool stale() const {
    std::error_code ec;
    const auto sz = std::filesystem::file_size(path, ec);
    if (ec || sz != size) {
      return true;
    }
    const auto mt = std::filesystem::last_write_time(path, ec);
    return ec || mt != mtime;
  }

  bool discoverNext() {
    if (!file.is_open()) {
      return false;
    }
    if (offsets.empty()) {
      file.clear();
      file.seekg(0);
    } else {
      file.clear();
      file.seekg(static_cast<std::streamoff>(offsets.back()));
      std::string skip;
      if (!std::getline(file, skip)) {
        return false;
      }
    }
    while (file) {
      const auto here = file.tellg();
      if (here < 0) {
        return false;
      }
      std::string line;
      if (!std::getline(file, line)) {
        return false;
      }
      if (isLammpsTimestep(line)) {
        offsets.push_back(static_cast<std::uint64_t>(here));
        return true;
      }
    }
    return false;
  }

  bool positionAt(int frame) {
    if (frame < 1 || !file.is_open()) {
      return false;
    }
    if (lastFrame + 1 == frame && lastFrame > 0) {
      const auto here = file.tellg();
      if (here >= 0 && static_cast<int>(offsets.size()) < frame) {
        offsets.push_back(static_cast<std::uint64_t>(here));
      }
      std::string line;
      if (!std::getline(file, line) || !isLammpsTimestep(line)) {
        return false;
      }
      lastFrame = frame;
      return true;
    }
    while (static_cast<int>(offsets.size()) < frame) {
      if (!discoverNext()) {
        return false;
      }
    }
    file.clear();
    file.seekg(static_cast<std::streamoff>(offsets[frame - 1]));
    std::string line;
    if (!std::getline(file, line) || !isLammpsTimestep(line)) {
      return false;
    }
    lastFrame = frame;
    return true;
  }

  int nframes() {
    while (discoverNext()) {
    }
    return static_cast<int>(offsets.size());
  }

  bool ensureIndexed(int frame) {
    if (frame < 1) {
      return false;
    }
    while (static_cast<int>(offsets.size()) < frame) {
      if (!discoverNext()) {
        return false;
      }
    }
    return true;
  }
};

std::mutex gDumpMu;
std::unordered_map<std::string, std::shared_ptr<LammpsDumpSession>> gDumps;

std::shared_ptr<LammpsDumpSession> sessionFor(const std::string &filename) {
  if (!gen::file_exists(filename)) {
    return nullptr;
  }
  std::lock_guard<std::mutex> lock(gDumpMu);
  auto &slot = gDumps[filename];
  if (!slot || slot->stale() || !slot->file.is_open()) {
    slot = std::make_shared<LammpsDumpSession>();
    if (!slot->open(filename)) {
      slot.reset();
      return nullptr;
    }
  }
  return slot;
}

void noteCoordColumn(int col, const std::string &tok, int &x, int &y,
                     int &z) {
  // LAMMPS ReaderNative (src/reader_native.cpp): take x if present,
  // otherwise the first of xs / xu / xsu in column order.
  if (tok == "x") {
    x = col;
  } else if ((tok == "xu" || tok == "xs" || tok == "xsu") && x < 0) {
    x = col;
  } else if (tok == "y") {
    y = col;
  } else if ((tok == "yu" || tok == "ys" || tok == "ysu") && y < 0) {
    y = col;
  } else if (tok == "z") {
    z = col;
  } else if ((tok == "zu" || tok == "zs" || tok == "zsu") && z < 0) {
    z = col;
  }
}

} // namespace

int sinp::nLammpsFrames(const std::string &filename) {
  auto sess = sessionFor(filename);
  if (sess == nullptr) {
    return 0;
  }
  std::lock_guard<std::mutex> lock(sess->mu);
  return sess->nframes();
}

void sinp::dropLammpsDumpIndex(const std::string &filename) {
  std::lock_guard<std::mutex> lock(gDumpMu);
  gDumps.erase(filename);
}

void sinp::forEachLammpsFrame(
    const std::string &filename, int first, int last, int typeFilter,
    const std::function<void(
        int, molSys::PointCloud<molSys::Point<double>, double> &)> &fn,
    int nThreads) {
  const int nframes = nLammpsFrames(filename);
  if (nframes <= 0) {
    return;
  }
  if (first < 1) {
    first = 1;
  }
  if (last <= 0 || last > nframes) {
    last = nframes;
  }
  if (first > last) {
    return;
  }

#ifdef SEAMS_HAS_OPENMP
  const int threads = nThreads > 0 ? nThreads : omp_get_max_threads();
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)             \
    if (threads > 1 && last > first)
#endif
  for (int frame = first; frame <= last; ++frame) {
    molSys::PointCloud<molSys::Point<double>, double> cloud;
    if (typeFilter > 0) {
      cloud = readLammpsTrjO(filename, frame, cloud, typeFilter);
    } else {
      cloud = readLammpsTrj(filename, frame, cloud);
    }
    fn(frame, cloud);
  }
}

/**
 * @details Get all the ring information, from the R.I.N.G.S. file. Each line
 * contains the IDs of the atoms in the ring. This is saved inside a vector of
 * vectors. Rings which have more than three consecutive water molecules are
 * discarded.
 */
std::vector<std::vector<int>> sinp::readBonds(std::string filename) {
  std::unique_ptr<std::ifstream> inpFile;
  inpFile = std::make_unique<std::ifstream>(filename);
  std::vector<std::vector<int>> bonds;
  std::string line;                // Current line being read in
  std::vector<std::string> tokens; // Vector containing word tokens
  std::vector<int> id;             // Vector for the IDs in the ring

  if (!(gen::file_exists(filename))) {
    std::cout << "Fatal Error: The bond file does not exist or you gave the "
                 "wrong path.\n";
    // Throw exception?
    return bonds;
  }

  // Format of the file:
  // 481 Bonds
  // 272    214    906   1361    388      1
  // 388   1361   1042   1548    237      1
  // 272   1536   1582   1701   1905      1

  if (inpFile->is_open()) {
    // ----------------------------------------------------------
    // At this point we know that the file is open
    std::getline((*inpFile), line); // Read in bonds
    // Run this until EOF or you reach the next timestep
    while (std::getline((*inpFile), line)) {
      // Read in lines and tokenize them into std::string words
      tokens = gen::tokenizer(line);

      id.clear();

      for (int i = 0; i < tokens.size(); i++) {
        id.push_back(std::stoi(tokens[i]));
      } // end of for, assigning values to id

      bonds.push_back(id);

    } // end of while, till EOF
    // ----------------------------------------------------------
  } // End of if file open statement

  inpFile->close();

  return bonds;
}
/**
 * @details Function for reading in an XYZ file
 */
molSys::PointCloud<molSys::Point<double>, double> sinp::readXYZ(std::string filename) {
  std::unique_ptr<std::ifstream> xyzFile;
  xyzFile = std::make_unique<std::ifstream>(filename);
  molSys::PointCloud<molSys::Point<double>, double>
      yCloud;
  std::string line;                // Current line being read in
  std::vector<std::string> tokens; // Vector containing word tokens
  std::vector<double> numbers;     // Vector containing type double numbers
  int nop = -1;                    // Number of atoms in targetFrame
  int iatom = 0;
  molSys::Point<double> iPoint; // Current point being read in from the file
  double xLo = 0.0, xHi = 0.0, yLo = 0.0, yHi = 0.0, zLo = 0.0, zHi = 0.0;

  if (!(gen::file_exists(filename))) {
    std::cout
        << "Fatal Error: The file does not exist or you gave the wrong path.\n";
    // Throw exception?
        throw std::runtime_error("Wrong filepath");
  }

  // --------
  // Before filling up the PointCloud, if the vectors are filled
  // empty them
  //yCloud = molSys::clearPointCloud(yCloud);
  // --------

  // Format of an XYZ file:
  //  1970
  // generated by VMD
  //  O         43.603500       16.926201       15.215700
  //  O         39.912601       14.775100       19.379200
  if (xyzFile->is_open()) {

    // ----------------------------------------------------------
    // At this point we know that the XYZ file is open

    // The first line contains the number of atoms
    std::getline((*xyzFile), line);
    nop = std::stoi(line);
    // Skip the comment line
    std::getline((*xyzFile), line);
    // Reserve memory for atom coordinates using nop read in
    yCloud.pts.reserve(nop);
    yCloud.nop = nop;
    // Run this until EOF or you reach the next timestep
    while (std::getline((*xyzFile), line)) {

      // Read in lines and tokenize them into std::string words and <double>
      // numbers
      tokens = gen::tokenizer(line);
      numbers = gen::tokenizerDouble(line);

      // Skip whitespace
      if (tokens.size() == 0) {
        continue;
      }

      // Put logic for checking atom type here later
      iPoint.type = 1; // Oxygen type; hard-coded!
      iPoint.x = std::stod(tokens[1]);
      iPoint.y = std::stod(tokens[2]);
      iPoint.z = std::stod(tokens[3]);
      if (yCloud.pts.empty()) {
        xLo = xHi = iPoint.x;
        yLo = yHi = iPoint.y;
        zLo = zHi = iPoint.z;
      } else {
        xLo = std::min(xLo, iPoint.x);
        xHi = std::max(xHi, iPoint.x);
        yLo = std::min(yLo, iPoint.y);
        yHi = std::max(yHi, iPoint.y);
        zLo = std::min(zLo, iPoint.z);
        zHi = std::max(zHi, iPoint.z);
      }
      iatom++;
      iPoint.molID = iatom;
      iPoint.atomID = iatom;
      yCloud.pts.push_back(iPoint);
      mapAtomIdToIndex(yCloud);
    } // end of while, looping through lines till EOF
    // ----------------------------------------------------------
  } // End of if file open statement

  xyzFile->close();

  if (yCloud.pts.size() == 1) {
    xHi = xLo + 10;
    yHi = yLo + 10;
    zHi = zLo + 10;
  } // for a single point in the system (never happens)

  yCloud.box = {xHi - xLo, yHi - yLo, zHi - zLo};
  yCloud.boxLow = {xLo, yLo, zLo};

  if (yCloud.pts.size() != yCloud.nop) {
    std::cout << "Atoms didn't get filled in properly.\n";
  }

  return yCloud;
}

// External Libraries

namespace {
enum class LammpsKeep { All, Type, TypeInSlice };

// Parse one snapshot. The session has already consumed ITEM: TIMESTEP.
// Stop at the next ITEM: TIMESTEP and seek back so the live cursor sits
// on that line, matching LAMMPS ReaderNative::read_time.
void parseLammpsFrameBody(
    std::ifstream &file,
    molSys::PointCloud<molSys::Point<double>, double> &yCloud, int typeFilter,
    LammpsKeep keep, bool isSlice, std::array<double, 3> coordLow,
    std::array<double, 3> coordHigh) {
  std::string line;
  std::vector<std::string> tokens;
  std::vector<double> numbers;
  std::vector<double> tilt;
  int nop = -1;
  bool readNOP = false;
  bool readBox = false;
  bool readAtoms = false;
  int xIndex = -1;
  int yIndex = -1;
  int zIndex = -1;
  int typeIndex = -1;
  int molIndex = 0;
  int atomIndex = 0;
  bool isTriclinic = false;
  int nKept = 0;
  molSys::Point<double> iPoint;

  while (true) {
    const auto mark = file.tellg();
    if (!std::getline(file, line)) {
      break;
    }
    if (isLammpsTimestep(line)) {
      if (mark >= 0) {
        file.clear();
        file.seekg(mark);
      }
      break;
    }

    tokens = gen::tokenizer(line);
    numbers = gen::tokenizerDouble(line);

    if (readNOP) {
      nop = std::stoi(line.data());
      readNOP = false;
      if (keep == LammpsKeep::All) {
        yCloud.pts.reserve(static_cast<std::size_t>(nop));
        yCloud.nop = nop;
      }
    }
    if (readBox) {
      if (!tokens.empty() && tokens[0] == "ITEM:") {
        readBox = false;
        if (isTriclinic) {
          for (std::size_t k = 0; k < tilt.size(); k++) {
            yCloud.box.push_back(tilt[k]);
          }
        }
      } else if (numbers.size() >= 2) {
        yCloud.box.push_back(numbers[1] - numbers[0]);
        yCloud.boxLow.push_back(numbers[0]);
        if (numbers.size() == 3) {
          isTriclinic = true;
          tilt.push_back(numbers[2]);
        }
      }
    }
    if (readAtoms && typeIndex >= 0 && xIndex >= 0 && yIndex >= 0 &&
        zIndex >= 0 &&
        numbers.size() >
            static_cast<std::size_t>(
                std::max({typeIndex, xIndex, yIndex, zIndex, molIndex,
                          atomIndex}))) {
      iPoint.type = static_cast<int>(numbers[typeIndex]);
      iPoint.molID = static_cast<int>(numbers[molIndex]);
      iPoint.atomID = static_cast<int>(numbers[atomIndex]);
      iPoint.x = numbers[xIndex];
      iPoint.y = numbers[yIndex];
      iPoint.z = numbers[zIndex];
      if (isSlice) {
        iPoint.inSlice = sinp::atomInSlice(iPoint.x, iPoint.y, iPoint.z,
                                           coordLow, coordHigh);
        if (keep == LammpsKeep::TypeInSlice && !iPoint.inSlice) {
          continue;
        }
      }
      if (keep != LammpsKeep::All && iPoint.type != typeFilter) {
        continue;
      }
      nKept++;
      yCloud.pts.push_back(iPoint);
      mapAtomIdToIndex(yCloud);
    }

    if (!tokens.empty() && tokens[0] == "ITEM:" && tokens.size() > 1) {
      if (tokens[1] == "NUMBER") {
        readNOP = true;
      } else if (tokens[1] == "BOX") {
        readBox = true;
      } else if (tokens[1] == "ATOMS") {
        readAtoms = true;
        xIndex = yIndex = zIndex = typeIndex = -1;
        molIndex = 0;
        atomIndex = 0;
        for (int i = 2; i < static_cast<int>(tokens.size()); i++) {
          if (tokens[i] == "type") {
            typeIndex = i - 2;
          } else if (tokens[i] == "mol") {
            molIndex = i - 2;
          } else if (tokens[i] == "id") {
            atomIndex = i - 2;
          } else {
            noteCoordColumn(i - 2, tokens[i], xIndex, yIndex, zIndex);
          }
        }
        if (molIndex == 0) {
          molIndex = atomIndex;
        }
      }
    }
  }

  if (keep != LammpsKeep::All) {
    yCloud.nop = static_cast<int>(yCloud.pts.size());
    if (yCloud.pts.size() != static_cast<std::size_t>(nKept)) {
      std::cout << "Atoms didn't get filled in properly.\n";
    }
  } else if (nop >= 0 && yCloud.pts.size() != static_cast<std::size_t>(yCloud.nop)) {
    std::cout << "Atoms didn't get filled in properly.\n";
  }
}

bool openCloneAt(const std::string &filename, int frame, std::ifstream &in) {
  auto sess = sessionFor(filename);
  if (sess == nullptr) {
    return false;
  }
  std::uint64_t off = 0;
  {
    std::lock_guard<std::mutex> lock(sess->mu);
    if (!sess->ensureIndexed(frame)) {
      return false;
    }
    off = sess->offsets[static_cast<std::size_t>(frame - 1)];
  }
  in.open(filename, std::ios::in | std::ios::binary);
  if (!in.is_open()) {
    return false;
  }
  in.seekg(static_cast<std::streamoff>(off));
  std::string line;
  return static_cast<bool>(std::getline(in, line)) && isLammpsTimestep(line);
}

void loadLammpsFrame(
    const std::string &filename, int targetFrame,
    molSys::PointCloud<molSys::Point<double>, double> &yCloud, int typeFilter,
    LammpsKeep keep, bool isSlice, std::array<double, 3> coordLow,
    std::array<double, 3> coordHigh) {
  yCloud = molSys::clearPointCloud(yCloud);
  if (!gen::file_exists(filename)) {
    std::cout
        << "Fatal Error: The file does not exist or you gave the wrong path.\n";
    yCloud.currentFrame = targetFrame;
    return;
  }
  auto sess = sessionFor(filename);
  if (sess == nullptr) {
    std::cout
        << "Fatal Error: The file does not exist or you gave the wrong path.\n";
    yCloud.currentFrame = targetFrame;
    return;
  }
  std::unique_lock<std::mutex> lock(sess->mu, std::try_to_lock);
  if (lock.owns_lock()) {
    if (!sess->positionAt(targetFrame)) {
      std::cout << "You entered a frame that doesn't exist.\n";
      yCloud.currentFrame = targetFrame;
      return;
    }
    parseLammpsFrameBody(sess->file, yCloud, typeFilter, keep, isSlice,
                         coordLow, coordHigh);
    yCloud.currentFrame = targetFrame;
    return;
  }
  std::ifstream clone;
  if (!openCloneAt(filename, targetFrame, clone)) {
    std::cout << "You entered a frame that doesn't exist.\n";
    yCloud.currentFrame = targetFrame;
    return;
  }
  parseLammpsFrameBody(clone, yCloud, typeFilter, keep, isSlice, coordLow,
                       coordHigh);
  yCloud.currentFrame = targetFrame;
}
} // namespace

/**
 * @details  Function for reading in a lammps file. Reads in a specified frame
 *  (frame number and not timestep value).
 * @param[in] filename The name of the lammps trajectory file to be read in
 * @param[in] targetFrame The frame number whose information will be read in
 * @param[out] yCloud The outputted PointCloud
 * @param[in] isSlice This decides whether a slice will be created or not
 * @param[in] coordLow Contains the lower limits of the slice, if a slice is to
 *  be created
 * @param[in] coordHigh Contains the upper limits of the slice, if a slice is to
 *  be created
 */
molSys::PointCloud<molSys::Point<double>, double>
sinp::readLammpsTrj(std::string filename, int targetFrame,
                    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                    bool isSlice, std::array<double, 3> coordLow,
                    std::array<double, 3> coordHigh) {
  loadLammpsFrame(filename, targetFrame, yCloud, -1, LammpsKeep::All, isSlice,
                  coordLow, coordHigh);
  return yCloud;
}

/**
 * @details Function for reading in a lammps file; and saves only the Oxygen
 * atoms. This is an overloaded function. The Oxygen atom ID must be specified.
 * @param[in] filename The name of the lammps trajectory file to be read in
 * @param[in] targetFrame The frame number whose information will be read in
 * @param[out] yCloud The outputted PointCloud
 * @param[in] typeO The type ID of the Oxygen atoms
 * @param[in] isSlice This decides whether a slice will be created or not
 * @param[in] coordLow Contains the lower limits of the slice, if a slice is to
 *  be created
 * @param[in] coordHigh Contains the upper limits of the slice, if a slice is
 *  to be created
 */
molSys::PointCloud<molSys::Point<double>, double>
sinp::readLammpsTrjO(std::string filename, int targetFrame,
                     molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                     int typeO, bool isSlice, std::array<double, 3> coordLow,
                     std::array<double, 3> coordHigh) {
  loadLammpsFrame(filename, targetFrame, yCloud, typeO, LammpsKeep::Type,
                  isSlice, coordLow, coordHigh);
  return yCloud;
}

/**
 * @details Function for reading in a lammps file; and saves only the atoms of
 * the desired type. Atoms which are not inside the slice or not of type I are
 * not saved at all This is an overloaded function. The type atom ID must be
 *  specified.
 * @param[in] filename The name of the lammps trajectory file to be read in
 * @param[in] targetFrame The frame number whose information will be read in
 * @param[out] yCloud The outputted PointCloud
 * @param[in] typeI The type ID of the desired type of atoms
 * @param[in] isSlice This decides whether a slice will be created or not
 * @param[in] coordLow Contains the lower limits of the slice, if a slice is to
 *  be created
 *  @param[in] coordHigh Contains the upper limits of the slice, if a slice is
 *  to be created
 */
molSys::PointCloud<molSys::Point<double>, double> sinp::readLammpsTrjreduced(
    std::string filename, int targetFrame,
    molSys::PointCloud<molSys::Point<double>, double> &yCloud, int typeI,
    bool isSlice, std::array<double, 3> coordLow,
    std::array<double, 3> coordHigh) {
  loadLammpsFrame(filename, targetFrame, yCloud, typeI,
                  LammpsKeep::TypeInSlice, isSlice, coordLow, coordHigh);
  return yCloud;
}

// ============================================================================
// Optional format readers (compiled only when dependencies are available)
// ============================================================================

#ifdef SEAMS_HAS_CHEMFILES
#include <chemfiles.hpp>

molSys::PointCloud<molSys::Point<double>, double>
sinp::readChemfiles(std::string filename, int targetFrame,
                    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                    int typeFilter) {
  try {
  chemfiles::Trajectory trajectory(filename);

  if (targetFrame < 1 ||
      static_cast<size_t>(targetFrame) > trajectory.nsteps()) {
    std::cerr << "Frame " << targetFrame << " does not exist in " << filename
              << " (has " << trajectory.nsteps() << " frames).\n";
    return yCloud;
  }

  yCloud = molSys::clearPointCloud(yCloud);

  auto frame = trajectory.read_step(static_cast<size_t>(targetFrame - 1));
  auto positions = frame.positions();
  auto cell = frame.cell();

  // Box dimensions
  yCloud.box = {cell.lengths()[0], cell.lengths()[1], cell.lengths()[2]};
  yCloud.boxLow = {0.0, 0.0, 0.0};

  auto &topology = frame.topology();

  for (size_t i = 0; i < frame.size(); i++) {
    // LAMMPS dumps carry numeric type identifiers, which chemfiles stores as
    // the atom type string with no atomic number; chemical formats carry
    // element names with an atomic number. Prefer a numeric type string,
    // fall back to the atomic number.
    const auto numeric = [](const std::string &s) {
      return !s.empty() &&
             s.find_first_not_of("0123456789") == std::string::npos;
    };
    int atomType = 0;
    if (numeric(topology[i].type())) {
      atomType = std::stoi(topology[i].type());
    } else if (numeric(topology[i].name())) {
      atomType = std::stoi(topology[i].name());
    } else {
      atomType = static_cast<int>(topology[i].atomic_number().value_or(1));
    }

    // Apply type filter if requested (-1 means accept all)
    if (typeFilter >= 0 && atomType != typeFilter) {
      continue;
    }

    molSys::Point<double> pt;
    pt.type = atomType;
    // chemfiles orders LAMMPS dump atoms by their id, so index + 1 recovers
    // the 1-based dump id for the contiguous case; other formats have no ID
    // notion beyond position
    pt.atomID = static_cast<int>(i) + 1;
    // Molecule IDs surface as chemfiles residues where the format has them
    const auto residue = topology.residue_for_atom(i);
    pt.molID = (residue && residue->id())
                   ? static_cast<int>(*residue->id())
                   : pt.atomID;
    pt.x = positions[i][0];
    pt.y = positions[i][1];
    pt.z = positions[i][2];

    yCloud.pts.push_back(pt);
    mapAtomIdToIndex(yCloud);
  }

  yCloud.nop = static_cast<int>(yCloud.pts.size());
  yCloud.currentFrame = targetFrame;

  return yCloud;
  } catch (const std::exception &e) {
    // chemfiles throws on unreadable or malformed files; report and return
    // the empty cloud instead of terminating the caller
    std::cerr << "chemfiles cannot read " << filename << ": " << e.what()
              << "\n";
    yCloud = molSys::clearPointCloud(yCloud);
    return yCloud;
  }
}
#endif // SEAMS_HAS_CHEMFILES

#ifdef SEAMS_HAS_READCON
// The C API rather than readcon-core.hpp: the C++ wrapper is a convenience
// layer over exactly these calls, and the C ABI is the surface cargo-c
// versions
#include <readcon-core.h>

molSys::PointCloud<molSys::Point<double>, double>
sinp::readCon(std::string filename, int targetFrame,
              molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  yCloud = molSys::clearPointCloud(yCloud);

  readcon::CConFrameIterator *frames =
      readcon::read_con_file_iterator(filename.c_str());
  if (frames == nullptr) {
    std::cerr << "Cannot open .con file " << filename << "\n";
    return yCloud;
  }

  int frameIdx = 0;
  while (readcon::RKRConFrame *handle =
             readcon::con_frame_iterator_next(frames)) {
    frameIdx++;
    if (frameIdx != targetFrame) {
      readcon::free_rkr_frame(handle);
      continue;
    }

    // Found the target frame; extract the transparent atom records
    readcon::CFrame *frame = readcon::rkr_frame_to_c_frame(handle);
    readcon::free_rkr_frame(handle);
    if (frame == nullptr) {
      std::cerr << "Cannot extract frame " << targetFrame << " from "
                << filename << "\n";
      break;
    }

    yCloud.box = {frame->cell[0], frame->cell[1], frame->cell[2]};
    yCloud.boxLow = {0.0, 0.0, 0.0};
    yCloud.pts.reserve(frame->num_atoms);

    for (size_t i = 0; i < frame->num_atoms; i++) {
      const readcon::CAtom &atom = frame->atoms[i];
      molSys::Point<double> pt;
      pt.type = static_cast<int>(atom.atomic_number);
      pt.atomID = static_cast<int>(atom.atom_id);
      pt.molID = pt.atomID;
      pt.x = atom.x;
      pt.y = atom.y;
      pt.z = atom.z;

      yCloud.pts.push_back(pt);
      mapAtomIdToIndex(yCloud);
    }
    readcon::free_c_frame(frame);

    yCloud.nop = static_cast<int>(yCloud.pts.size());
    yCloud.currentFrame = targetFrame;
    readcon::free_con_frame_iterator(frames);
    return yCloud;
  }

  if (yCloud.pts.empty()) {
    std::cerr << "Frame " << targetFrame << " not found in " << filename
              << " (has " << frameIdx << " frames).\n";
  }
  readcon::free_con_frame_iterator(frames);
  return yCloud;
}
#endif // SEAMS_HAS_READCON
