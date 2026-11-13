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

#include <algorithm>
#include <atomic>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <vector>
#include <filesystem>
#include <fstream>
#include <generic.hpp>
#include <map>
#include <sstream>
#include <neighbours.hpp>
#include <memory>
#include <mutex>
#include <seams_input.hpp>
#include <string_view>
#include <thread>
#include <sys/stat.h>
#include <unistd.h>
#include <ext/stdio_filebuf.h>
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

bool isLammpsItem(const std::string &line) {
  return line.size() >= 5 && line.compare(0, 5, "ITEM:") == 0;
}

// Locale-independent field scan. istringstream >> double copies the line
// and walks a facet; from_chars is the C++17 path Lemire measured at
// gigabyte-per-second scale on the same decimal grammar LAMMPS writes.
constexpr int kMaxDumpFields = 32;

int parseDumpFields(std::string_view line, double *out, int cap) {
  const char *p = line.data();
  const char *const end = p + line.size();
  int n = 0;
  while (p < end && n < cap) {
    while (p < end && (*p == ' ' || *p == '\t' || *p == '\r')) {
      ++p;
    }
    if (p >= end) {
      break;
    }
    double value = 0;
    // from_chars(double) is the fast C-locale scan. libc++ on the macOS 11
    // conda SDK deletes it, and that build keeps strtod.
#if defined(__cpp_lib_to_chars) && __cpp_lib_to_chars >= 201611L
    const auto parsed = std::from_chars(p, end, value);
    if (parsed.ec != std::errc{}) {
      while (p < end && *p != ' ' && *p != '\t') {
        ++p;
      }
      continue;
    }
    out[n++] = value;
    p = parsed.ptr;
#else
    char buf[64];
    const auto ntok = std::min<std::size_t>(static_cast<std::size_t>(end - p),
                                            sizeof(buf) - 1);
    std::memcpy(buf, p, ntok);
    buf[ntok] = '\0';
    char *endp = nullptr;
    value = std::strtod(buf, &endp);
    if (endp == buf) {
      while (p < end && *p != ' ' && *p != '\t') {
        ++p;
      }
      continue;
    }
    out[n++] = value;
    p += (endp - buf);
#endif
  }
  return n;
}

bool parseDumpInt(std::string_view line, int &out) {
  const char *p = line.data();
  const char *const end = p + line.size();
  while (p < end && (*p == ' ' || *p == '\t')) ++p;
  const auto parsed = std::from_chars(p, end, out);
  return parsed.ec == std::errc{};
}

// One POSIX descriptor. The stream reads that descriptor and nothing else.
class OwnedFdStream {
public:
  explicit OwnedFdStream(int fd)
      : fd_(fd), buf_(fd, std::ios::in | std::ios::binary), in_(&buf_) {}
  std::istream &stream() { return in_; }
  int descriptor() const { return fd_; }
  bool failed() const { return in_.bad(); }

private:
  int fd_;
  __gnu_cxx::stdio_filebuf<char> buf_;
  std::istream in_;
};

bool streamReadFailed(std::istream &in) { return in.bad(); }

// Live dump position, matching LAMMPS ReaderNative (rerun keeps FILE* at
// the next snapshot) and chemfiles LAMMPSTrajectory::read_next. Random
// access seeks a cached ITEM: TIMESTEP offset table.
struct LammpsDumpSession {
  std::mutex mu;
  std::string path;
  std::unique_ptr<OwnedFdStream> io;
  std::vector<std::uint64_t> offsets;
  int lastFrame = 0;
  bool fullyIndexed = false;
  std::filesystem::file_time_type mtime{};
  std::uintmax_t size{0};

  int descriptor() const { return io ? io->descriptor() : -1; }
  std::istream &stream() { return io->stream(); }

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
    const int raw = ::open(filename.c_str(), O_RDONLY);
    if (raw < 0) {
      return false;
    }
    io = std::make_unique<OwnedFdStream>(raw);
    offsets.clear();
    lastFrame = 0;
    fullyIndexed = false;
    return static_cast<bool>(stream());
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

  // A private copy, then a memchr walk of ITEM: TIMESTEP. An mmap of
  // the file faults with SIGBUS when another writer truncates or
  // replaces it before the pages are read. pread stops at the bytes
  // that still belong to this descriptor.
  bool mmapIndexAll() {
    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) {
      return false;
    }
    struct stat st {};
    if (::fstat(fd, &st) != 0 || st.st_size <= 0) {
      ::close(fd);
      return false;
    }
    const auto nbytes = static_cast<std::size_t>(st.st_size);
    std::vector<char> bytes(nbytes);
    std::size_t got = 0;
    while (got < nbytes) {
      const ssize_t n = ::pread(fd, bytes.data() + got, nbytes - got,
                               static_cast<off_t>(got));
      if (n < 0) {
        if (errno == EINTR) {
          continue;
        }
        ::close(fd);
        return false;
      }
      if (n == 0) {
        break;
      }
      got += static_cast<std::size_t>(n);
    }
    ::close(fd);
    offsets.clear();
    const char *const begin = bytes.data();
    const char *const end = begin + got;
    constexpr char kNeedle[] = "ITEM: TIMESTEP";
    constexpr std::size_t kN = 14;
    const char *p = begin;
    while (p + kN <= end) {
      const void *hit = std::memchr(p, 'I', static_cast<std::size_t>(end - p));
      if (hit == nullptr) {
        break;
      }
      p = static_cast<const char *>(hit);
      if (p + kN <= end && std::memcmp(p, kNeedle, kN) == 0 &&
          (p == begin || p[-1] == '\n')) {
        offsets.push_back(static_cast<std::uint64_t>(p - begin));
        p += kN;
      } else {
        ++p;
      }
    }
    fullyIndexed = true;
    return !offsets.empty();
  }

  bool discoverNext() {
    if (offsets.empty()) {
      stream().clear();
      stream().seekg(0);
    } else {
      stream().clear();
      stream().seekg(static_cast<std::streamoff>(offsets.back()));
      std::string skip;
      if (!std::getline(stream(), skip)) {
        return false;
      }
    }
    while (stream()) {
      const auto here = stream().tellg();
      if (here < 0) {
        return false;
      }
      std::string line;
      if (!std::getline(stream(), line)) {
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
    if (frame < 1) {
      return false;
    }
    if (lastFrame + 1 == frame && lastFrame > 0) {
      const auto here = stream().tellg();
      if (here >= 0 && static_cast<int>(offsets.size()) < frame) {
        offsets.push_back(static_cast<std::uint64_t>(here));
      }
      std::string line;
      if (!std::getline(stream(), line) || !isLammpsTimestep(line)) {
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
    stream().clear();
    stream().seekg(static_cast<std::streamoff>(offsets[static_cast<std::size_t>(frame - 1)]));
    std::string line;
    if (!std::getline(stream(), line) || !isLammpsTimestep(line)) {
      return false;
    }
    lastFrame = frame;
    return true;
  }

  int nframes() {
    if (!fullyIndexed) {
      mmapIndexAll();
    }
    if (!fullyIndexed) {
      while (discoverNext()) {
      }
      fullyIndexed = true;
    }
    return static_cast<int>(offsets.size());
  }

  bool ensureIndexed(int frame) {
    if (frame < 1) {
      return false;
    }
    if (!fullyIndexed && static_cast<int>(offsets.size()) < frame) {
      mmapIndexAll();
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
  if (!slot || slot->stale()) {
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
  try {
    auto sess = sessionFor(filename);
    if (sess == nullptr) {
      return 0;
    }
    std::lock_guard<std::mutex> lock(sess->mu);
    return sess->nframes();
  } catch (const std::ios_base::failure &) {
    return 0;
  }
}

void sinp::dropLammpsDumpIndex(const std::string &filename) {
  std::lock_guard<std::mutex> lock(gDumpMu);
  gDumps.erase(filename);
}

void sinp::forEachLammpsFrame(
    const std::string &filename, int first, int last, int typeFilter,
    const std::function<void(
        int, molSys::PointCloud<molSys::Point<double>, double> &)> &fn,
    int nThreads, int part, int parts) {
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
  parts = std::max(parts, 1);
  if (part < 0 || part >= parts) {
    return;
  }
  std::vector<int> frames;
  for (int round = 0; first + round * parts <= last; ++round) {
    const int slot = round % 2 == 0 ? part : parts - 1 - part;
    const int frame = first + round * parts + slot;
    if (frame <= last) {
      frames.push_back(frame);
    }
  }
  const int count = static_cast<int>(frames.size());

  std::mutex errorMu;
  std::atomic<bool> failed{false};
  std::string error;
#ifdef SEAMS_HAS_OPENMP
  const int threads = nThreads > 0 ? nThreads : omp_get_max_threads();
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)             \
    if (threads > 1 && count > 1)
#endif
  for (int k = 0; k < count; ++k) {
    if (failed.load()) {
      continue;
    }
    const int frame = frames[static_cast<std::size_t>(k)];
    try {
      molSys::PointCloud<molSys::Point<double>, double> cloud;
      if (typeFilter > 0) {
        cloud = readLammpsTrjO(filename, frame, cloud, typeFilter);
      } else {
        cloud = readLammpsTrj(filename, frame, cloud);
      }
      fn(frame, cloud);
    } catch (const std::exception &ex) {
      std::lock_guard<std::mutex> lock(errorMu);
      if (error.empty()) {
        error = ex.what();
      }
      failed.store(true);
    }
  }
  if (!error.empty()) {
    throw std::runtime_error(error);
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
      bool ok = true;
      for (int i = 0; i < static_cast<int>(tokens.size()); i++) {
        try {
          id.push_back(std::stoi(tokens[static_cast<std::size_t>(i)]));
        } catch (...) {
          ok = false;
          break;
        }
      }
      if (!ok || id.empty()) {
        continue;
      }
      bonds.push_back(id);

    } // end of while, till EOF
    // The comment above this function shows a RINGS line whose last
    // column is a constant 0/1 flag, not an atom ID. Drop that column
    // only when every row carries the same flag and is wider than a pair.
    bool dropFlag = !bonds.empty();
    int flag = 0;
    for (std::size_t row = 0; row < bonds.size(); row++) {
      const auto &ids = bonds[row];
      if (ids.size() < 3 || (ids.back() != 0 && ids.back() != 1)) {
        dropFlag = false;
        break;
      }
      if (row == 0) {
        flag = ids.back();
      } else if (ids.back() != flag) {
        dropFlag = false;
        break;
      }
    }
    if (dropFlag) {
      for (auto &ids : bonds) {
        ids.pop_back();
      }
    }
    // ----------------------------------------------------------
  } // End of if file open statement

  inpFile->close();

  return bonds;
}
namespace {

// Extended XYZ stores the cell as Lattice="ax ay az bx by bz cx cy cz".
// A reduced cell becomes the LAMMPS box (lengths, then xy xz yz). The
// coordinate span is not a cell.
bool latticeCell(const std::string &comment, std::vector<double> &box,
                 std::vector<double> &boxLow) {
  const auto key = comment.find("Lattice=");
  if (key == std::string::npos) {
    return false;
  }
  const auto open = comment.find_first_of("\"'", key + 8);
  if (open == std::string::npos) {
    throw std::runtime_error("XYZ Lattice is not a cell.");
  }
  const char quote = comment[open];
  const auto close = comment.find(quote, open + 1);
  if (close == std::string::npos) {
    throw std::runtime_error("XYZ Lattice is not a cell.");
  }
  std::istringstream in(comment.substr(open + 1, close - open - 1));
  double cell[9];
  for (double &value : cell) {
    if (!(in >> value)) {
      throw std::runtime_error("XYZ Lattice is not a cell.");
    }
  }
  const double ax = cell[0];
  const double ay = cell[1];
  const double az = cell[2];
  const double bx = cell[3];
  const double by = cell[4];
  const double bz = cell[5];
  const double cx = cell[6];
  const double cy = cell[7];
  const double cz = cell[8];
  const double scale = std::max(std::abs(ax), std::max(std::abs(by), std::abs(cz)));
  const double tol = 1e-8 * (scale > 1.0 ? scale : 1.0);
  if (!(ax > 0.0 && by > 0.0 && cz > 0.0) || std::abs(ay) > tol ||
      std::abs(az) > tol || std::abs(bz) > tol) {
    throw std::runtime_error("XYZ Lattice is not a reduced cell.");
  }
  // Same reduced frame as setBoxFromLammpsEdges: a along x, b in the
  // xy plane. The stored box is the bound span plus xy, xz, yz, and
  // the origin shifts when a tilt is negative.
  const double lx = ax;
  const double ly = by;
  const double lz = cz;
  const double xy = bx;
  const double xz = cx;
  const double yz = cy;
  const double xmin = std::min(std::min(0.0, xy), std::min(xz, xy + xz));
  const double xmax = std::max(std::max(0.0, xy), std::max(xz, xy + xz));
  const double ymin = std::min(0.0, yz);
  const double ymax = std::max(0.0, yz);
  const bool tilt =
      std::abs(xy) > tol || std::abs(xz) > tol || std::abs(yz) > tol;
  if (!tilt) {
    box = {lx, ly, lz};
    boxLow = {0.0, 0.0, 0.0};
    return true;
  }
  box = {lx + xmax - xmin, ly + ymax - ymin, lz, xy, xz, yz};
  boxLow = {xmin, ymin, 0.0};
  return true;
}

} // namespace

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
  bool haveLattice = false;

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
    if (!std::getline((*xyzFile), line)) {
      xyzFile->close();
      throw std::runtime_error("XYZ header is not an atom count.");
    }
    try {
      nop = std::stoi(line);
    } catch (...) {
      xyzFile->close();
      throw std::runtime_error("XYZ header is not an atom count.");
    }
    if (nop < 0) {
      xyzFile->close();
      throw std::runtime_error("XYZ header is not an atom count.");
    }
    // The comment may carry an extended-XYZ Lattice.
    std::string comment;
    std::getline((*xyzFile), comment);
    try {
      haveLattice = latticeCell(comment, yCloud.box, yCloud.boxLow);
    } catch (...) {
      xyzFile->close();
      throw;
    }
    // Reserve memory for atom coordinates using nop read in
    yCloud.pts.reserve(static_cast<std::size_t>(nop));
    yCloud.nop = nop;
    // Run this until EOF or you reach the next timestep
    std::map<std::string, int> elementType;
    // Stop at the declared count. A later frame is the next record,
    // not extra atoms of this one.
    while (static_cast<int>(yCloud.pts.size()) < nop &&
           std::getline((*xyzFile), line)) {

      // Read in lines and tokenize them into std::string words and <double>
      // numbers
      tokens = gen::tokenizer(line);
      numbers = gen::tokenizerDouble(line);

      // Skip whitespace
      if (tokens.size() == 0) {
        continue;
      }

      const std::size_t coordinateOffset = tokens.size() >= 4 ? 1 : 0;
      if (tokens.size() < coordinateOffset + 3) {
        continue;
      }

      // One element stays type 1. A second element takes the next type,
      // in order of first appearance, so Si then N is 1 then 2.
      if (coordinateOffset == 1) {
        const std::string &symbol = tokens[0];
        char *end = nullptr;
        std::strtod(symbol.c_str(), &end);
        const bool numeric = end != symbol.c_str() && *end == '\0';
        if (numeric) {
          iPoint.type = 1;
        } else {
          auto found = elementType.find(symbol);
          if (found == elementType.end()) {
            const int id = static_cast<int>(elementType.size()) + 1;
            found = elementType.emplace(symbol, id).first;
          }
          iPoint.type = found->second;
        }
      } else {
        iPoint.type = 1;
      }
      try {
        iPoint.x = std::stod(tokens[coordinateOffset]);
        iPoint.y = std::stod(tokens[coordinateOffset + 1]);
        iPoint.z = std::stod(tokens[coordinateOffset + 2]);
      } catch (...) {
        xyzFile->close();
        throw std::runtime_error("XYZ coordinate is not a number.");
      }
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

  if (!haveLattice && yCloud.pts.size() == 1) {
    xHi = xLo + 10;
    yHi = yLo + 10;
    zHi = zLo + 10;
  } // for a single point in the system (never happens)

  if (!haveLattice) {
    yCloud.box = {xHi - xLo, yHi - yLo, zHi - zLo};
    yCloud.boxLow = {xLo, yLo, zLo};
  }
  yCloud.currentFrame = 1;

  if (yCloud.pts.size() != static_cast<std::size_t>(yCloud.nop)) {
    // A short or over-long frame is not a structure. Callers that kept
    // the header count used to index past pts, and callers that rewrote
    // nop used to report the fragment as one frame.
    throw std::runtime_error(
        "XYZ frame is truncated: expected " + std::to_string(yCloud.nop) +
        " atoms, read " + std::to_string(yCloud.pts.size()) + ".");
  }

  return yCloud;
}

// External Libraries

namespace {
enum class LammpsKeep { All, Type, TypeInSlice };

// The descriptor was opened with the stream, so it names the same inode.
// fstat's result is checked. A failure is a change, and the stat is not used.
bool dumpChangedDuringRead(int fd, const std::string &path,
                           const struct stat &opened) {
  struct stat now {};
  if (::fstat(fd, &now) != 0) {
    return true;
  }
  if (now.st_size != opened.st_size) {
    return true;
  }
  struct stat atPath {};
  if (stat(path.c_str(), &atPath) != 0) {
    return true;
  }
  return atPath.st_ino != now.st_ino;
}

// A writer can append in the gap between the last byte and this check.
// One short poll keeps that snapshot. A file that stays the same size
// and inode through the poll is finished, and a short frame is an error.
bool dumpStillGrowing(int fd, const std::string &path, const struct stat &opened) {
  if (dumpChangedDuringRead(fd, path, opened)) {
    return true;
  }
  const auto deadline =
      std::chrono::steady_clock::now() + std::chrono::milliseconds(50);
  while (std::chrono::steady_clock::now() < deadline) {
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
    if (dumpChangedDuringRead(fd, path, opened)) {
      return true;
    }
  }
  return false;
}

// Parse one snapshot. The session has already consumed ITEM: TIMESTEP.
// Stop at the next ITEM: TIMESTEP and seek back so the live position sits
// on that line, matching LAMMPS ReaderNative::read_time.
void parseLammpsFrameBody(
    std::istream &file, int fd, const std::string &path,
    molSys::PointCloud<molSys::Point<double>, double> &yCloud, int typeFilter,
    LammpsKeep keep, bool isSlice, std::array<double, 3> coordLow,
    std::array<double, 3> coordHigh) {
  struct stat opened {};
  const bool haveOpened = fstat(fd, &opened) == 0;
  std::string line;
  std::vector<std::string> tokens;
  std::vector<double> tilt;
  int nop = -1;
  bool readNOP = false;
  bool readBox = false;
  bool readAtoms = false;
  int xIndex = -1;
  int yIndex = -1;
  int zIndex = -1;
  int typeIndex = -1;
  int molIndex = -1;
  int atomIndex = -1;
  bool sawId = false;
  bool sawMol = false;
  int nextSyntheticId = 1;
  bool isTriclinic = false;
  int nLines = 0;
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

    const bool itemLine = isLammpsItem(line);
    if (itemLine) {
      tokens = gen::tokenizer(line);
    } else {
      tokens.clear();
    }

    if (readNOP && !itemLine) {
      if (parseDumpInt(line, nop) && nop >= 0) {
        yCloud.pts.reserve(static_cast<std::size_t>(nop));
        yCloud.idIndexMap.reserve(static_cast<std::size_t>(nop));
        if (keep == LammpsKeep::All) {
          yCloud.nop = nop;
        }
      }
      readNOP = false;
    }
    if (readBox) {
      if (itemLine) {
        readBox = false;
        if (isTriclinic) {
          for (std::size_t k = 0; k < tilt.size(); k++) {
            yCloud.box.push_back(tilt[k]);
          }
        }
      } else {
        double fields[kMaxDumpFields];
        const int n = parseDumpFields(line, fields, kMaxDumpFields);
        if (n >= 2) {
          yCloud.box.push_back(fields[1] - fields[0]);
          yCloud.boxLow.push_back(fields[0]);
          if (n >= 3) {
            isTriclinic = true;
            tilt.push_back(fields[2]);
          }
        }
      }
    }
    if (readAtoms && !itemLine && typeIndex >= 0 && xIndex >= 0 &&
        yIndex >= 0 && zIndex >= 0) {
      if (line.find_first_not_of(" \t\r") != std::string::npos) {
        nLines += 1;
      }
      double fields[kMaxDumpFields];
      const int n = parseDumpFields(line, fields, kMaxDumpFields);
      int need = std::max({typeIndex, xIndex, yIndex, zIndex});
      if (sawId) {
        need = std::max(need, atomIndex);
      }
      if (sawMol) {
        need = std::max(need, molIndex);
      }
      // n is the field count, and need is a 0-based column, so the line
      // must hold at least need + 1 fields.
      if (n > need) {
        iPoint.type = static_cast<int>(fields[typeIndex]);
        // A dump with no id column cannot reuse the type column: every
        // oxygen would share one ID and the last atom would win the map.
        if (sawId) {
          iPoint.atomID = static_cast<int>(fields[atomIndex]);
        } else {
          iPoint.atomID = nextSyntheticId++;
        }
        // mol in column 0 is a real column. A missing mol column used to
        // be detected by molIndex staying 0, which overwrote that column.
        if (sawMol) {
          iPoint.molID = static_cast<int>(fields[molIndex]);
        } else {
          iPoint.molID = iPoint.atomID;
        }
        iPoint.x = fields[xIndex];
        iPoint.y = fields[yIndex];
        iPoint.z = fields[zIndex];
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
        yCloud.pts.push_back(iPoint);
        mapAtomIdToIndex(yCloud);
      }
    }

    if (itemLine && tokens.size() > 1) {
      if (tokens[1] == "NUMBER") {
        readNOP = true;
      } else if (tokens[1] == "BOX") {
        readBox = true;
      } else if (tokens[1] == "ATOMS") {
        readAtoms = true;
        xIndex = yIndex = zIndex = typeIndex = -1;
        molIndex = -1;
        atomIndex = -1;
        sawId = false;
        sawMol = false;
        nextSyntheticId = 1;
        for (int i = 2; i < static_cast<int>(tokens.size()); i++) {
          if (tokens[static_cast<std::size_t>(i)] == "type") {
            typeIndex = i - 2;
          } else if (tokens[static_cast<std::size_t>(i)] == "mol") {
            molIndex = i - 2;
            sawMol = true;
          } else if (tokens[static_cast<std::size_t>(i)] == "id") {
            atomIndex = i - 2;
            sawId = true;
          } else {
            noteCoordColumn(i - 2, tokens[static_cast<std::size_t>(i)], xIndex, yIndex, zIndex);
          }
        }
      }
    }
  }

  // Fewer atom lines than the declared count is not a frame when the
  // file is finished. A count that does not fit in int leaves nop unset
  // and is not this case. A dump still being written changes under the
  // read; that snapshot keeps the atoms it has.
  // A read error is the stream state (errno on the owned descriptor),
  // not a finished short file. A stable size and inode is a finished file.
  if (nop >= 0 && nLines < nop && !streamReadFailed(file) &&
      !(haveOpened && dumpStillGrowing(fd, path, opened))) {
    throw std::runtime_error(
        "LAMMPS dump frame is truncated: expected " + std::to_string(nop) +
        " atoms, read " + std::to_string(nLines) + ".");
  }

  // A frame cut short (a dump still being written) keeps the atoms it has;
  // the warning goes to stderr, which every rank of a frame split keeps
  if (keep != LammpsKeep::All) {
    yCloud.nop = static_cast<int>(yCloud.pts.size());
  } else if (yCloud.pts.size() !=
             static_cast<std::size_t>(yCloud.nop < 0 ? 0 : yCloud.nop)) {
    // Reaching this else-if means the stored count and the atom vector
    // disagree. A missing count leaves the stored nop at 0, so an empty
    // vector never enters, and the warning condition is always true here.
    std::cerr << "Atoms didn't get filled in properly.\n";
    yCloud.nop = static_cast<int>(yCloud.pts.size());
  }
}

struct CloneRead {
  std::unique_ptr<OwnedFdStream> io;
};

bool openCloneAt(const std::string &filename, int frame, CloneRead &out) {
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
  const int raw = ::open(filename.c_str(), O_RDONLY);
  if (raw < 0) {
    return false;
  }
  out.io = std::make_unique<OwnedFdStream>(raw);
  out.io->stream().seekg(static_cast<std::streamoff>(off));
  std::string line;
  if (!std::getline(out.io->stream(), line) || !isLammpsTimestep(line)) {
    out.io.reset();
    return false;
  }
  return true;
}

bool parseLooseDouble(const char *&p, const char *end, double &out) {
  while (p < end && (*p == ' ' || *p == '\t')) {
    ++p;
  }
  if (p >= end) {
    return false;
  }
#if defined(__cpp_lib_to_chars) && __cpp_lib_to_chars >= 201611L
  const auto parsed = std::from_chars(p, end, out);
  if (parsed.ec != std::errc{}) {
    return false;
  }
  p = parsed.ptr;
  return true;
#else
  char buf[64];
  const auto ntok =
      std::min<std::size_t>(static_cast<std::size_t>(end - p), sizeof(buf) - 1);
  std::memcpy(buf, p, ntok);
  buf[ntok] = '\0';
  char *endp = nullptr;
  out = std::strtod(buf, &endp);
  if (endp == buf) {
    return false;
  }
  p += endp - buf;
  return true;
#endif
}

bool parseLooseInt(const char *&p, const char *end, int &out) {
  while (p < end && (*p == ' ' || *p == '\t')) {
    ++p;
  }
  if (p >= end) {
    return false;
  }
  const auto parsed = std::from_chars(p, end, out);
  if (parsed.ec != std::errc{}) {
    return false;
  }
  p = parsed.ptr;
  return true;
}

// id type x y z on one orthorhombic frame. Anything else returns false and
// the streaming parser reads the frame. The map is filled in one pass so a
// 1e5-atom dump does not rehash and does not cross the iostream per line.
bool parseMappedIdTypeXyz(const std::string &filename, std::uint64_t frameAt,
                          molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                          int typeFilter, LammpsKeep keep, std::uint64_t &nextAt) {
  if (keep == LammpsKeep::TypeInSlice) {
    return false;
  }
  const int fd = ::open(filename.c_str(), O_RDONLY);
  if (fd < 0) {
    return false;
  }
  struct stat st {};
  if (::stat(filename.c_str(), &st) != 0 || st.st_size <= 0) {
    ::close(fd);
    return false;
  }
  const auto nbytes = static_cast<std::size_t>(st.st_size);
  if (frameAt >= nbytes) {
    ::close(fd);
    return false;
  }
  // pread, not mmap: a writer that truncates the dump faults an mmap with SIGBUS.
  std::vector<char> bytes(nbytes);
  std::size_t got = 0;
  while (got < nbytes) {
    const ssize_t nread = ::pread(fd, bytes.data() + got, nbytes - got,
                                  static_cast<off_t>(got));
    if (nread < 0) {
      if (errno == EINTR) {
        continue;
      }
      ::close(fd);
      return false;
    }
    if (nread == 0) {
      break;
    }
    got += static_cast<std::size_t>(nread);
  }
  ::close(fd);
  if (got <= frameAt) {
    return false;
  }
  bytes.resize(got);
  const char *const begin = bytes.data();
  const char *const end = begin + bytes.size();
  const char *p = begin + frameAt;
  auto nextLine = [&](const char *&s, const char *&e) {
    s = p;
    while (p < end && *p != '\n') {
      ++p;
    }
    e = p;
    if (e > s && e[-1] == '\r') {
      --e;
    }
    if (p < end && *p == '\n') {
      ++p;
    }
  };
  const char *s = nullptr;
  const char *e = nullptr;
  auto startsWith = [](const char *s, const char *e, const char *lit) {
    const std::size_t n = std::strlen(lit);
    return static_cast<std::size_t>(e - s) >= n && std::memcmp(s, lit, n) == 0;
  };
  nextLine(s, e);
  if (!startsWith(s, e, "ITEM: TIMESTEP")) {
    return false;
  }
  nextLine(s, e);
  nextLine(s, e);
  if (!startsWith(s, e, "ITEM: NUMBER")) {
    return false;
  }
  nextLine(s, e);
  int nop = 0;
  const char *np = s;
  if (!parseLooseInt(np, e, nop) || nop < 0) {
    return false;
  }
  nextLine(s, e);
  if (!startsWith(s, e, "ITEM: BOX") ||
      std::string_view(s, static_cast<std::size_t>(e - s)).find("xy") !=
          std::string_view::npos) {
    return false;
  }
  double span[3];
  double low[3];
  for (int axis = 0; axis < 3; ++axis) {
    nextLine(s, e);
    const char *bp = s;
    double lo = 0;
    double hi = 0;
    if (!parseLooseDouble(bp, e, lo) || !parseLooseDouble(bp, e, hi)) {
      return false;
    }
    low[axis] = lo;
    span[axis] = hi - lo;
  }
  nextLine(s, e);
  while (e > s && (e[-1] == ' ' || e[-1] == '\t')) {
    --e;
  }
  if (std::string_view(s, static_cast<std::size_t>(e - s)) !=
      "ITEM: ATOMS id type x y z") {
    return false;
  }
  yCloud.pts.clear();
  yCloud.pts.reserve(static_cast<std::size_t>(nop));
  yCloud.idIndexMap.clear();
  yCloud.idIndexMap.reserve(static_cast<std::size_t>(nop));
  yCloud.box = {span[0], span[1], span[2]};
  yCloud.boxLow = {low[0], low[1], low[2]};
  const bool filter = keep != LammpsKeep::All;
  int kept = 0;
  bool failed = false;
  for (int i = 0; i < nop; ++i) {
    if (p >= end) {
      failed = true;
      break;
    }
    nextLine(s, e);
    if (startsWith(s, e, "ITEM:")) {
      failed = true;
      break;
    }
    const char *ap = s;
    int id = 0;
    int typ = 0;
    double x = 0;
    double y = 0;
    double z = 0;
    if (!parseLooseInt(ap, e, id) || !parseLooseInt(ap, e, typ) ||
        !parseLooseDouble(ap, e, x) || !parseLooseDouble(ap, e, y) ||
        !parseLooseDouble(ap, e, z)) {
      failed = true;
      break;
    }
    if (filter && typ != typeFilter) {
      continue;
    }
    molSys::Point<double> point;
    point.atomID = id;
    point.molID = id;
    point.type = typ;
    point.x = x;
    point.y = y;
    point.z = z;
    yCloud.pts.push_back(point);
    if (!yCloud.idIndexMap.emplace(id, kept).second) {
      std::cerr << "Warning: duplicate atom ID " << id
                << " in the input frame; the ID-to-index map keeps only the "
                   "last atom read with this ID\n";
      yCloud.idIndexMap[id] = kept;
    }
    ++kept;
  }
  if (failed) {
    // The snapshot was taken at st.st_size. An append or a replace during
    // the copy changes that size or the inode. The atoms already parsed
    // are the frame; falling through would open the grown file, see a
    // stable short count, and report a finished truncated dump.
    struct stat later {};
    const bool moved = ::stat(filename.c_str(), &later) != 0 ||
                       later.st_size != st.st_size || later.st_ino != st.st_ino;
    if (moved && kept > 0) {
      yCloud.nop = kept;
      nextAt = static_cast<std::uint64_t>(p - begin);
      return true;
    }
    yCloud.pts.clear();
    yCloud.idIndexMap.clear();
    yCloud.box.clear();
    yCloud.boxLow.clear();
    return false;
  }
  if (keep == LammpsKeep::All && kept != nop) {
    std::cerr << "Atoms didn't get filled in properly.\n";
  }
  yCloud.nop = kept;
  nextAt = static_cast<std::uint64_t>(p - begin);
  return true;
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
    std::uint64_t nextAt = 0;
    if (!isSlice && targetFrame >= 1 &&
        static_cast<std::size_t>(targetFrame) <= sess->offsets.size() &&
        parseMappedIdTypeXyz(
            filename, sess->offsets[static_cast<std::size_t>(targetFrame - 1)],
            yCloud, typeFilter, keep, nextAt)) {
      sess->stream().clear();
      sess->stream().seekg(static_cast<std::streamoff>(nextAt));
      yCloud.currentFrame = targetFrame;
      return;
    }
    parseLammpsFrameBody(sess->stream(), sess->descriptor(), filename, yCloud,
                         typeFilter, keep, isSlice, coordLow, coordHigh);
    yCloud.currentFrame = targetFrame;
    return;
  }
  CloneRead clone;
  if (!openCloneAt(filename, targetFrame, clone)) {
    std::cout << "You entered a frame that doesn't exist.\n";
    yCloud.currentFrame = targetFrame;
    return;
  }
  parseLammpsFrameBody(clone.io->stream(), clone.io->descriptor(), filename, yCloud,
                       typeFilter, keep, isSlice, coordLow, coordHigh);
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
 * @details One LAMMPS dump frame, one atom type. The type argument is
 *  any LAMMPS type (the O in the name is historical). isSlice sets
 *  inSlice and does not drop atoms. nop is the type-filtered count.
 * @param[in] filename The name of the lammps trajectory file to be read in
 * @param[in] targetFrame The frame number whose information will be read in
 * @param[out] yCloud The outputted PointCloud
 * @param[in] typeO The LAMMPS type ID to keep
 * @param[in] isSlice If true, set inSlice from the box; do not drop
 * @param[in] coordLow Lower limits of the slice. An axis with lo == hi
 *  is unconstrained.
 * @param[in] coordHigh Upper limits of the slice
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
 * @details One LAMMPS dump frame, one atom type, dropping atoms
 *  outside the slice when isSlice is true. nop is the kept count.
 * @param[in] filename The name of the lammps trajectory file to be read in
 * @param[in] targetFrame The frame number whose information will be read in
 * @param[out] yCloud The outputted PointCloud
 * @param[in] typeI The LAMMPS type ID to keep
 * @param[in] isSlice If true, drop atoms outside the box
 * @param[in] coordLow Lower limits of the slice. An axis with lo == hi
 *  is unconstrained.
 * @param[in] coordHigh Upper limits of the slice
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

#ifdef SEAMS_HAS_READCON
#include <readcon-core.hpp>
#include <numbers>
#ifdef SEAMS_HAS_READCON_DB
#include <readcon-db.h>
#endif

namespace {

bool conAnglesOrtho(const double angles[3]) {
  for (int k = 0; k < 3; k++) {
    if (std::abs(angles[k] - 90.0) > 1e-8) {
      return false;
    }
  }
  return true;
}

// CON lengths and angles are the ASE / eOn crystallographic basis.
// An orthorhombic cell stays three bound spans so the pair loop does
// not take the triclinic path. A shear is stored as dump spans plus
// xy, xz, yz, the same box a tilt dump already carries. The conversion
// does not need libminimage.
void applyConCell(molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                  const double cell[3], const double angles[3]) {
  const double aLen = cell[0];
  const double bLen = cell[1];
  const double cLen = cell[2];
  yCloud.boxLow = {0.0, 0.0, 0.0};
  if (!(aLen > 0.0 && bLen > 0.0 && cLen > 0.0)) {
    yCloud.box = {aLen, bLen, cLen};
    return;
  }
  constexpr double deg = std::numbers::pi / 180.0;
  const double cosA = std::cos(angles[0] * deg);
  const double cosB = std::cos(angles[1] * deg);
  const double cosG = std::cos(angles[2] * deg);
  const double sinG = std::sin(angles[2] * deg);
  if (std::abs(sinG) < 1e-18) {
    yCloud.box = {aLen, bLen, cLen};
    return;
  }
  const double xy = bLen * cosG;
  const double ly = bLen * sinG;
  const double xz = cLen * cosB;
  const double yz = cLen * (cosA - cosB * cosG) / sinG;
  const double cz2 = cLen * cLen - xz * xz - yz * yz;
  if (!(cz2 > 0.0)) {
    yCloud.box = {aLen, bLen, cLen};
    return;
  }
  const double lz = std::sqrt(cz2);
  if (conAnglesOrtho(angles)) {
    yCloud.box = {aLen, ly, lz};
    return;
  }
  const double xmin = std::min(std::min(0.0, xy), std::min(xz, xy + xz));
  const double xmax = std::max(std::max(0.0, xy), std::max(xz, xy + xz));
  const double ymin = std::min(0.0, yz);
  const double ymax = std::max(0.0, yz);
  yCloud.box = {aLen + xmax - xmin, ly + ymax - ymin, lz, xy, xz, yz};
  yCloud.boxLow = {xmin, ymin, 0.0};
}

bool cloudFromHandle(
    const readcon::RKRConFrame *handle, int frameNo,
    molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  readcon::CFrame *frame = readcon::rkr_frame_to_c_frame(handle);
  if (frame == nullptr) {
    return false;
  }
  applyConCell(yCloud, frame->cell, frame->angles);
  uintptr_t nxyz = 0;
  const double *xyz = readcon::rkr_frame_xyz_f64(handle, &nxyz);
  // n is the atom count and the buffer is row-major xyz. Matching only
  // y and z of atom 0 accepts a column-major buffer when that atom has
  // y == z and the next x happens to equal them.
  const bool rowMajor =
      xyz != nullptr && nxyz == frame->num_atoms &&
      (frame->num_atoms == 0 ||
       (xyz[0] == frame->atoms[0].x && xyz[1] == frame->atoms[0].y &&
        xyz[2] == frame->atoms[0].z &&
        (frame->num_atoms < 2 || xyz[3] == frame->atoms[1].x)));
  yCloud.pts.reserve(frame->num_atoms);
  for (size_t i = 0; i < frame->num_atoms; i++) {
    const readcon::CAtom &atom = frame->atoms[i];
    molSys::Point<double> pt;
    pt.type = static_cast<int>(atom.atomic_number);
    pt.atomID = static_cast<int>(atom.atom_id);
    pt.molID = pt.atomID;
    if (rowMajor) {
      pt.x = xyz[3 * i];
      pt.y = xyz[3 * i + 1];
      pt.z = xyz[3 * i + 2];
    } else {
      pt.x = atom.x;
      pt.y = atom.y;
      pt.z = atom.z;
    }
    yCloud.pts.push_back(pt);
    mapAtomIdToIndex(yCloud);
  }
  readcon::free_c_frame(frame);
  yCloud.nop = static_cast<int>(yCloud.pts.size());
  yCloud.currentFrame = frameNo;
  return true;
}

#ifdef SEAMS_HAS_READCON_DB
// shards.json is the manifest readcon-db's campaign commands write
// (shard-init, drain, join). The C ABI opens one LMDB env, so a root
// is fanned out by that layout instead of a second parser.
uint32_t campaignShardCount(const std::filesystem::path &manifest) {
  std::ifstream in(manifest);
  if (!in) {
    return 0;
  }
  std::string text((std::istreambuf_iterator<char>(in)),
                   std::istreambuf_iterator<char>());
  const auto key = text.find("\"n_shards\"");
  if (key == std::string::npos) {
    return 0;
  }
  const auto colon = text.find(':', key);
  if (colon == std::string::npos) {
    return 0;
  }
  const unsigned long n = std::strtoul(text.c_str() + colon + 1, nullptr, 10);
  if (n == 0 || n > 4096) {
    return 0;
  }
  return static_cast<uint32_t>(n);
}

std::string campaignDest(const std::string &corpus, std::uint64_t trajId) {
  const auto manifest = std::filesystem::path(corpus) / "shards.json";
  const uint32_t n = campaignShardCount(manifest);
  if (n == 0) {
    return corpus;
  }
  char name[32];
  std::snprintf(name, sizeof name, "shard_%04u",
                static_cast<unsigned>(trajId % n));
  return (std::filesystem::path(corpus) / name).string();
}

std::vector<std::string> campaignEnvs(const std::string &corpus) {
  const auto manifest = std::filesystem::path(corpus) / "shards.json";
  const uint32_t n = campaignShardCount(manifest);
  std::vector<std::string> roots;
  if (n == 0) {
    roots.push_back(corpus);
    return roots;
  }
  for (uint32_t i = 0; i < n; i++) {
    char name[32];
    std::snprintf(name, sizeof name, "shard_%04u", i);
    const auto dir = std::filesystem::path(corpus) / name;
    if (std::filesystem::is_regular_file(dir / "data.mdb")) {
      roots.push_back(dir.string());
    }
  }
  return roots;
}
#endif

} // namespace

molSys::PointCloud<molSys::Point<double>, double>
sinp::readCon(std::string filename, int targetFrame,
              molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  yCloud = molSys::clearPointCloud(yCloud);
  if (targetFrame < 1) {
    return yCloud;
  }
  try {
    readcon::ConFrameIterator frames(filename);
    std::optional<readcon::ConFrame> frame =
        frames.nth(static_cast<std::size_t>(targetFrame - 1));
    if (!frame) {
      std::cerr << "Frame " << targetFrame << " not found in " << filename
                << "\n";
      return yCloud;
    }
    if (!cloudFromHandle(frame->get_handle(), targetFrame, yCloud)) {
      std::cerr << "Cannot extract frame " << targetFrame << " from "
                << filename << "\n";
      yCloud = molSys::clearPointCloud(yCloud);
    }
  } catch (const std::exception &e) {
    std::cerr << "Cannot open .con file " << filename << ": " << e.what()
              << "\n";
    yCloud = molSys::clearPointCloud(yCloud);
  }
  return yCloud;
}

#ifdef SEAMS_HAS_READCON_DB
int sinp::ingestConCorpus(const std::string &corpus, std::uint64_t trajId,
                          const std::string &conPath) {
  const std::string dest = campaignDest(corpus, trajId);
  std::error_code ec;
  std::filesystem::create_directories(dest, ec);
  if (ec) {
    std::cerr << "Cannot create corpus " << dest << ": " << ec.message()
              << "\n";
    return -1;
  }
  size_t id = 0;
  if (rkrdb_open(dest.c_str(), &id) != RKRDB_OK) {
    std::cerr << "Cannot open corpus " << dest << "\n";
    return -1;
  }
  uint32_t nframes = 0;
  const int st =
      rkrdb_append_trajectory(id, trajId, conPath.c_str(), &nframes);
  rkrdb_close(id);
  if (st != RKRDB_OK) {
    std::cerr << "Cannot ingest " << conPath << " into " << dest << "\n";
    return -1;
  }
  return static_cast<int>(nframes);
}

molSys::PointCloud<molSys::Point<double>, double> sinp::readConCorpus(
    const std::string &corpus, const std::string &symbol,
    const std::string &formula, int targetFrame,
    molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  yCloud = molSys::clearPointCloud(yCloud);
  if (targetFrame < 1) {
    return yCloud;
  }
  struct Hit {
    std::string root;
    std::uint64_t traj = 0;
    std::uint32_t frame = 0;
  };
  std::vector<Hit> hits;
  const char *sym = symbol.empty() ? nullptr : symbol.c_str();
  const char *form = formula.empty() ? nullptr : formula.c_str();
  for (const std::string &root : campaignEnvs(corpus)) {
    size_t id = 0;
    if (rkrdb_open_readonly(root.c_str(), &id) != RKRDB_OK) {
      continue;
    }
    const int st = rkrdb_select_campaign(
        id, -1, sym, 0, UINT32_MAX, form, 0.0, 0.0, 0, 0.0, 0.0, 0, nullptr,
        0, 0, 0, 0);
    if (st == RKRDB_OK) {
      const int n = rkrdb_result_count(id);
      for (int i = 0; i < n; i++) {
        Hit hit;
        hit.root = root;
        if (rkrdb_result_key(id, static_cast<size_t>(i), &hit.traj,
                             &hit.frame) == RKRDB_OK) {
          hits.push_back(hit);
        }
      }
    }
    rkrdb_close(id);
  }
  std::sort(hits.begin(), hits.end(), [](const Hit &a, const Hit &b) {
    if (a.traj != b.traj) {
      return a.traj < b.traj;
    }
    return a.frame < b.frame;
  });
  if (targetFrame > static_cast<int>(hits.size())) {
    std::cerr << "Frame " << targetFrame << " not found in corpus " << corpus
              << " (" << hits.size() << " hits).\n";
    return yCloud;
  }
  const Hit &hit = hits[static_cast<std::size_t>(targetFrame - 1)];
  size_t id = 0;
  if (rkrdb_open_readonly(hit.root.c_str(), &id) != RKRDB_OK) {
    std::cerr << "Cannot reopen corpus " << hit.root << "\n";
    return yCloud;
  }
  std::string text(1u << 20, '\0');
  int n = rkrdb_get_frame_text(id, hit.traj, hit.frame, text.data(),
                               text.size());
  if (n == RKRDB_ERR) {
    text.assign(8u << 20, '\0');
    n = rkrdb_get_frame_text(id, hit.traj, hit.frame, text.data(), text.size());
  }
  rkrdb_close(id);
  if (n < 0) {
    std::cerr << "Cannot read frame text from " << hit.root << "\n";
    return yCloud;
  }
  text.resize(static_cast<std::size_t>(n));
  readcon::CConFrameIterator *frames =
      readcon::read_con_string_iterator(text.c_str());
  if (frames == nullptr) {
    std::cerr << "readcon-core rejected the corpus blob\n";
    return yCloud;
  }
  readcon::RKRConFrame *handle = readcon::con_frame_iterator_next(frames);
  if (handle == nullptr || !cloudFromHandle(handle, targetFrame, yCloud)) {
    std::cerr << "Cannot decode corpus frame " << targetFrame << "\n";
    yCloud = molSys::clearPointCloud(yCloud);
  }
  if (handle != nullptr) {
    readcon::free_rkr_frame(handle);
  }
  readcon::free_con_frame_iterator(frames);
  return yCloud;
}
#endif // SEAMS_HAS_READCON_DB
#endif // SEAMS_HAS_READCON
