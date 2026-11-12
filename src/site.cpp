//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <site.hpp>

#include <franzblau.hpp>
#include <generic.hpp>
#include <mol_sys.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

namespace site {

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

std::string trim(std::string_view s) {
  const auto begin = s.find_first_not_of(" \t\n\r");
  if (begin == std::string_view::npos) {
    return {};
  }
  const auto end = s.find_last_not_of(" \t\n\r");
  return std::string(s.substr(begin, end - begin + 1));
}

Kind kindFromName(std::string_view name) {
  if (name == "unspecified") {
    return Kind::unspecified;
  }
  if (name == "cationHead") {
    return Kind::cationHead;
  }
  if (name == "anion") {
    return Kind::anion;
  }
  if (name == "tail") {
    return Kind::tail;
  }
  if (name == "donorH") {
    return Kind::donorH;
  }
  if (name == "acceptor") {
    return Kind::acceptor;
  }
  if (name == "polar") {
    return Kind::polar;
  }
  if (name == "apolar") {
    return Kind::apolar;
  }
  if (name == "waterO") {
    return Kind::waterO;
  }
  if (name == "waterH") {
    return Kind::waterH;
  }
  if (name == "solvent") {
    return Kind::solvent;
  }
  throw std::invalid_argument("unknown site kind '" + std::string(name) + "'");
}

Family familyFromName(std::string_view name) {
  if (name == "waterIce") {
    return Family::waterIce;
  }
  if (name == "ionicLiquid") {
    return Family::ionicLiquid;
  }
  if (name == "moltenSalt") {
    return Family::moltenSalt;
  }
  if (name == "des") {
    return Family::des;
  }
  if (name == "electrolyte") {
    return Family::electrolyte;
  }
  if (name == "confinedIL") {
    return Family::confinedIL;
  }
  if (name == "confinedWater") {
    return Family::confinedWater;
  }
  if (name == "networkFormer") {
    return Family::networkFormer;
  }
  throw std::invalid_argument("unknown site family '" + std::string(name) +
                              "'");
}

bool isIonKind(Kind k) {
  return k == Kind::cationHead || k == Kind::anion;
}

bool matchesKind(Kind have, Kind want) {
  if (want == Kind::polar) {
    return have == Kind::cationHead || have == Kind::anion ||
           have == Kind::polar;
  }
  if (want == Kind::apolar) {
    return have == Kind::tail || have == Kind::apolar;
  }
  return have == want;
}

int indexOfAtom(const Cloud &src, int atomID) {
  const auto it = src.idIndexMap.find(atomID);
  if (it != src.idIndexMap.end()) {
    return it->second;
  }
  const int n = static_cast<int>(src.pts.size());
  for (int i = 0; i < n; ++i) {
    if (src.pts[static_cast<std::size_t>(i)].atomID == atomID) {
      return i;
    }
  }
  return -1;
}

} // namespace

Kind Table::of(const molSys::Point<double> &p) const {
  if (const auto it = atomOverride.find(p.atomID); it != atomOverride.end()) {
    return it->second;
  }
  return ofType(p.type);
}

Kind Table::ofType(int typeId) const {
  if (const auto it = typeToKind.find(typeId); it != typeToKind.end()) {
    return it->second;
  }
  return Kind::unspecified;
}

std::vector<int> indicesOf(const Cloud &yCloud, const Table &table,
                           Kind kind) {
  std::vector<int> out;
  const int n = static_cast<int>(yCloud.pts.size());
  out.reserve(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    if (matchesKind(table.of(yCloud.pts[static_cast<std::size_t>(i)]), kind)) {
      out.push_back(i);
    }
  }
  return out;
}

int lammpsTypeOfKind(const Table &table, Kind kind) {
  int found = 0;
  int typeId = -1;
  for (const auto &[id, mapped] : table.typeToKind) {
    if (mapped == kind) {
      ++found;
      typeId = id;
    }
  }
  if (found != 1) {
    throw std::runtime_error("site kind does not map to a unique LAMMPS type");
  }
  return typeId;
}

Cloud ionCloud(const Cloud &src, const Table &table) {
  Cloud out;
  out.box = src.box;
  out.boxLow = src.boxLow;
  out.currentFrame = src.currentFrame;

  auto &mutableSrc = const_cast<Cloud &>(src);
  const auto molMap = molSys::createMolIDAtomIDMultiMap(mutableSrc);

  std::vector<int> molOrder;
  std::unordered_set<int> seenMol;
  molOrder.reserve(src.pts.size());
  for (const auto &pt : src.pts) {
    if (seenMol.insert(pt.molID).second) {
      molOrder.push_back(pt.molID);
    }
  }

  for (const int molID : molOrder) {
    std::vector<int> ions;
    const auto range = molMap.equal_range(molID);
    for (auto it = range.first; it != range.second; ++it) {
      const int idx = indexOfAtom(src, it->second);
      if (idx < 0) {
        continue;
      }
      if (isIonKind(table.of(src.pts[static_cast<std::size_t>(idx)]))) {
        ions.push_back(idx);
      }
    }
    if (ions.empty()) {
      continue;
    }
    std::sort(ions.begin(), ions.end());
    const Kind ionKind = table.of(src.pts[static_cast<std::size_t>(ions.front())]);
    ions.erase(std::remove_if(ions.begin(), ions.end(),
                              [&](int idx) {
                                return table.of(src.pts[static_cast<std::size_t>(
                                           idx)]) != ionKind;
                              }),
               ions.end());
    if (ions.empty()) {
      continue;
    }

    molSys::Point<double> vertex = src.pts[static_cast<std::size_t>(ions.front())];
    vertex.type = (ionKind == Kind::cationHead) ? 1 : 2;
    vertex.molID = molID;
    if (ions.size() > 1) {
      const int ref = ions.front();
      double sx = src.pts[static_cast<std::size_t>(ref)].x;
      double sy = src.pts[static_cast<std::size_t>(ref)].y;
      double sz = src.pts[static_cast<std::size_t>(ref)].z;
      for (std::size_t k = 1; k < ions.size(); ++k) {
        const auto dr = gen::relDist(src, ions[k], ref);
        sx += src.pts[static_cast<std::size_t>(ref)].x + dr[0];
        sy += src.pts[static_cast<std::size_t>(ref)].y + dr[1];
        sz += src.pts[static_cast<std::size_t>(ref)].z + dr[2];
      }
      const double inv = 1.0 / static_cast<double>(ions.size());
      vertex.x = sx * inv;
      vertex.y = sy * inv;
      vertex.z = sz * inv;
    }

    const int outIdx = static_cast<int>(out.pts.size());
    out.idIndexMap[vertex.atomID] = outIdx;
    out.pts.push_back(std::move(vertex));
  }

  out.nop = static_cast<int>(out.pts.size());
  return out;
}

Table parseSiteSpec(std::string_view spec) {
  Table table;
  std::size_t start = 0;
  while (start <= spec.size()) {
    const std::size_t comma = spec.find(',', start);
    const auto raw =
        spec.substr(start, comma == std::string_view::npos
                               ? std::string_view::npos
                               : comma - start);
    start = (comma == std::string_view::npos) ? spec.size() + 1 : comma + 1;
    const std::string token = trim(raw);
    if (token.empty()) {
      continue;
    }
    const auto eq = token.find('=');
    if (eq == std::string::npos) {
      throw std::invalid_argument("site spec token '" + token +
                                  "' is not key=value");
    }
    const std::string key = trim(token.substr(0, eq));
    const std::string val = trim(token.substr(eq + 1));
    if (key.empty() || val.empty()) {
      throw std::invalid_argument("site spec token '" + token +
                                  "' is not key=value");
    }
    if (key == "family") {
      table.family = familyFromName(val);
      continue;
    }
    std::size_t consumed = 0;
    int typeId = 0;
    try {
      typeId = std::stoi(key, &consumed);
    } catch (const std::exception &) {
      throw std::invalid_argument("site spec type '" + key +
                                  "' is not an integer");
    }
    if (consumed != key.size()) {
      throw std::invalid_argument("site spec type '" + key +
                                  "' is not an integer");
    }
    table.typeToKind[typeId] = kindFromName(val);
  }
  return table;
}

IonEnvironment
ionEnvironment(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
               const std::vector<bool> &iceFlag, const std::vector<int> &ionIndices,
               int waterType, double cutoff) {
  if (waterType == 0) {
    return ionEnvironment(yCloud, iceFlag, ionIndices, std::vector<int>{}, cutoff);
  }
  return ionEnvironment(yCloud, iceFlag, ionIndices, std::vector<int>{waterType}, cutoff);
}

IonEnvironment
ionEnvironment(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
               const std::vector<bool> &iceFlag, const std::vector<int> &ionIndices,
               const std::vector<int> &waterTypes, double cutoff) {
  IonEnvironment out;
  const int n = yCloud.nop;
  std::vector<char> isIon(static_cast<std::size_t>(std::max(n, 0)), 0);
  for (int i : ionIndices) {
    if (i >= 0 && i < n) {
      isIon[static_cast<std::size_t>(i)] = 1;
    }
  }
  const double cut2 = cutoff * cutoff;
  const gen::FracBox frame = gen::makeFracBox(yCloud);
  for (int i : ionIndices) {
    if (i < 0 || i >= n) {
      continue;
    }
    int shell = 0;
    int labelled = 0;
    std::vector<int> members;
    for (int j = 0; j < n; j++) {
      if (j == i || isIon[static_cast<std::size_t>(j)]) {
        continue;
      }
      if (!waterTypes.empty() &&
          std::find(waterTypes.begin(), waterTypes.end(),
                    yCloud.pts[static_cast<std::size_t>(j)].type) == waterTypes.end()) {
        continue;
      }
      if (gen::periodicDistSq(frame, yCloud, i, j) >= cut2) {
        continue;
      }
      ++shell;
      members.push_back(j);
      if (static_cast<std::size_t>(j) < iceFlag.size() && iceFlag[static_cast<std::size_t>(j)]) {
        ++labelled;
      }
    }
    IonState state = IonState::liquid;
    if (shell > 0 && labelled == shell) {
      state = IonState::ice;
      ++out.nIce;
    } else if (labelled > 0) {
      state = IonState::front;
      ++out.nFront;
    } else {
      ++out.nLiquid;
    }
    out.ion.push_back(i);
    out.shell.push_back(shell);
    out.iceFraction.push_back(shell > 0 ? static_cast<double>(labelled) / shell : 0.0);
    out.state.push_back(state);
    out.members.push_back(std::move(members));
  }
  return out;
}

std::vector<int> shellRingCensus(const std::vector<std::vector<int>> &rings,
                                 const std::vector<int> &shell, int maxRingSize) {
  std::vector<int> census(static_cast<std::size_t>(std::max(maxRingSize, 0)) + 1, 0);
  std::unordered_set<int> inShell(shell.begin(), shell.end());
  for (const auto &ring : rings) {
    if (ring.empty() || static_cast<int>(ring.size()) > maxRingSize) {
      continue;
    }
    for (int a : ring) {
      if (inShell.count(a)) {
        census[ring.size()] += 1;
        break;
      }
    }
  }
  return census;
}

std::vector<std::vector<int>> cappedRings(const std::vector<std::vector<int>> &rings,
                                          const std::vector<int> &shell,
                                          int maxRingSize) {
  std::vector<std::vector<int>> out;
  const std::unordered_set<int> inShell(shell.begin(), shell.end());
  for (const auto &ring : rings) {
    if (ring.empty() || static_cast<int>(ring.size()) > maxRingSize) {
      continue;
    }
    if (std::all_of(ring.begin(), ring.end(),
                    [&](int a) { return inShell.count(a) != 0; })) {
      out.push_back(ring);
    }
  }
  return out;
}

std::vector<std::vector<int>> brokenRings(const std::vector<std::vector<int>> &nList,
                                          const std::vector<int> &shell,
                                          int maxRingSize) {
  std::vector<std::vector<int>> out;
  const int n = static_cast<int>(nList.size());
  std::vector<int> members;
  for (int w : shell) {
    if (w >= 0 && w < n) {
      members.push_back(w);
    }
  }
  std::sort(members.begin(), members.end());
  members.erase(std::unique(members.begin(), members.end()), members.end());
  if (members.empty() || maxRingSize < 3) {
    return out;
  }
  const int ion = n;
  std::vector<std::vector<int>> aug;
  aug.reserve(static_cast<std::size_t>(n) + 1);
  aug.assign(nList.begin(), nList.end());
  std::vector<int> ionRow = {ion};
  for (int w : members) {
    auto &row = aug[static_cast<std::size_t>(w)];
    if (row.empty()) {
      row.push_back(w);
    }
    row.push_back(ion);
    ionRow.push_back(w);
  }
  aug.push_back(std::move(ionRow));

  // A ring through the ion has its lowest-indexed member within
  // maxRingSize / 2 hops of the ion, so only those sources are searched.
  std::vector<char> sources(aug.size(), 0);
  std::vector<int> dist(aug.size(), -1);
  std::vector<int> frontier = {ion};
  dist[static_cast<std::size_t>(ion)] = 0;
  sources[static_cast<std::size_t>(ion)] = 1;
  for (int depth = 1; depth <= maxRingSize / 2 && !frontier.empty(); depth++) {
    std::vector<int> next;
    for (int u : frontier) {
      const auto &row = aug[static_cast<std::size_t>(u)];
      for (std::size_t m = 1; m < row.size(); m++) {
        const int v = row[m];
        if (v >= 0 && v <= ion && dist[static_cast<std::size_t>(v)] < 0) {
          dist[static_cast<std::size_t>(v)] = depth;
          sources[static_cast<std::size_t>(v)] = 1;
          next.push_back(v);
        }
      }
    }
    frontier.swap(next);
  }
  for (auto &ring : primitive::ringNetwork(aug, maxRingSize, sources)) {
    if (std::find(ring.begin(), ring.end(), ion) != ring.end()) {
      out.push_back(std::move(ring));
    }
  }
  return out;
}


namespace {
// minimum-image displacement r - p for a point r and a reference p
std::array<double, 3>
minImage(const molSys::PointCloud<molSys::Point<double>, double> &yCloud, double rx,
         double ry, double rz, double px, double py, double pz) {
  if (yCloud.box.size() >= 6) {
    return gen::triclinicMinImage(yCloud, rx, ry, rz, px, py, pz);
  }
  std::array<double, 3> dr = {rx - px, ry - py, rz - pz};
  for (int k = 0; k < 3; k++) {
    const double L = yCloud.box[static_cast<std::size_t>(k)];
    if (L > 0.0) {
      dr[static_cast<std::size_t>(k)] -= L * std::round(dr[static_cast<std::size_t>(k)] / L);
    }
  }
  return dr;
}
} // namespace

std::array<double, 3>
periodicCentroid(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                 const std::vector<int> &atoms) {
  std::array<double, 3> c = {0.0, 0.0, 0.0};
  if (atoms.empty()) {
    return c;
  }
  const auto &p0 = yCloud.pts[static_cast<std::size_t>(atoms.front())];
  for (int a : atoms) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(a)];
    const auto dr = minImage(yCloud, p.x, p.y, p.z, p0.x, p0.y, p0.z);
    for (int k = 0; k < 3; k++) {
      c[static_cast<std::size_t>(k)] += dr[static_cast<std::size_t>(k)];
    }
  }
  const double inv = 1.0 / static_cast<double>(atoms.size());
  c[0] = p0.x + c[0] * inv;
  c[1] = p0.y + c[1] * inv;
  c[2] = p0.z + c[2] * inv;
  return c;
}

GuestOccupancy
guestOccupancy(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
               const std::vector<std::vector<int>> &cages,
               const std::vector<int> &guestIndices, double radius) {
  GuestOccupancy out;
  out.guestsPerCage.assign(cages.size(), 0);
  std::vector<std::array<double, 3>> centres;
  centres.reserve(cages.size());
  for (const auto &cage : cages) {
    centres.push_back(periodicCentroid(yCloud, cage));
  }
  const double r2max = radius * radius;
  for (int g : guestIndices) {
    if (g < 0 || g >= yCloud.nop) {
      continue;
    }
    const auto &p = yCloud.pts[static_cast<std::size_t>(g)];
    int best = -1;
    double bestSq = r2max;
    for (std::size_t c = 0; c < centres.size(); c++) {
      const auto dr = minImage(yCloud, p.x, p.y, p.z, centres[c][0], centres[c][1], centres[c][2]);
      const double d2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
      if (d2 <= bestSq) {
        bestSq = d2;
        best = static_cast<int>(c);
      }
    }
    out.cageOfGuest.push_back(best);
    out.centreDistance.push_back(best < 0 ? -1.0 : std::sqrt(bestSq));
    if (best < 0) {
      ++out.free;
    } else {
      ++out.guestsPerCage[static_cast<std::size_t>(best)];
    }
  }
  for (int n : out.guestsPerCage) {
    if (n > 0) {
      ++out.occupied;
    }
    if (n > 1) {
      ++out.multiply;
    }
  }
  int maxG = 0;
  for (int n : out.guestsPerCage) {
    maxG = std::max(maxG, n);
  }
  out.occupancyHistogram.assign(static_cast<std::size_t>(maxG) + 1, 0);
  for (int n : out.guestsPerCage) {
    ++out.occupancyHistogram[static_cast<std::size_t>(n)];
  }
  return out;
}

namespace {

using Vec3 = std::array<double, 3>;

Vec3 sub3(const Vec3 &a, const Vec3 &b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

Vec3 addScaled(const Vec3 &a, const Vec3 &d, double s) {
  return {a[0] + s * d[0], a[1] + s * d[1], a[2] + s * d[2]};
}

Vec3 crossP(const Vec3 &a, const Vec3 &b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

double dotP(const Vec3 &a, const Vec3 &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double distSq(const Vec3 &a, const Vec3 &b) {
  const auto d = sub3(a, b);
  return dotP(d, d);
}

// Squared distance between segments p1q1 and p2q2 (Ericson, Real-Time
// Collision Detection, 5.1.9).
double segmentDistSq(const Vec3 &p1, const Vec3 &q1, const Vec3 &p2,
                     const Vec3 &q2) {
  constexpr double tiny = 1e-300;
  const auto d1 = sub3(q1, p1);
  const auto d2 = sub3(q2, p2);
  const auto r = sub3(p1, p2);
  const double a = dotP(d1, d1);
  const double e = dotP(d2, d2);
  const double f = dotP(d2, r);
  double s = 0.0;
  double t = 0.0;
  if (a <= tiny && e <= tiny) {
    return dotP(r, r);
  }
  if (a <= tiny) {
    t = std::clamp(f / e, 0.0, 1.0);
  } else {
    const double c = dotP(d1, r);
    if (e <= tiny) {
      s = std::clamp(-c / a, 0.0, 1.0);
    } else {
      const double b = dotP(d1, d2);
      const double denom = a * e - b * b;
      s = denom > 0.0 ? std::clamp((b * f - c * e) / denom, 0.0, 1.0) : 0.0;
      t = (b * s + f) / e;
      if (t < 0.0) {
        t = 0.0;
        s = std::clamp(-c / a, 0.0, 1.0);
      } else if (t > 1.0) {
        t = 1.0;
        s = std::clamp((b - c) / a, 0.0, 1.0);
      }
    }
  }
  return distSq(addScaled(p1, d1, s), addScaled(p2, d2, t));
}

// Squared distance from p to triangle abc (Ericson, 5.1.5); a triangle of
// zero area is measured as its three edges.
double triangleDistSq(const Vec3 &p, const Vec3 &a, const Vec3 &b,
                      const Vec3 &c) {
  const auto ab = sub3(b, a);
  const auto ac = sub3(c, a);
  const auto n = crossP(ab, ac);
  const double scale = std::max(dotP(ab, ab), dotP(ac, ac));
  if (dotP(n, n) <= 1e-24 * scale * scale) {
    return std::min({segmentDistSq(p, p, a, b), segmentDistSq(p, p, b, c),
                     segmentDistSq(p, p, c, a)});
  }
  const auto ap = sub3(p, a);
  const double d1 = dotP(ab, ap);
  const double d2 = dotP(ac, ap);
  if (d1 <= 0.0 && d2 <= 0.0) {
    return distSq(p, a);
  }
  const auto bp = sub3(p, b);
  const double d3 = dotP(ab, bp);
  const double d4 = dotP(ac, bp);
  if (d3 >= 0.0 && d4 <= d3) {
    return distSq(p, b);
  }
  const double vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    return distSq(p, addScaled(a, ab, d1 / (d1 - d3)));
  }
  const auto cp = sub3(p, c);
  const double d5 = dotP(ab, cp);
  const double d6 = dotP(ac, cp);
  if (d6 >= 0.0 && d5 <= d6) {
    return distSq(p, c);
  }
  const double vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    return distSq(p, addScaled(a, ac, d2 / (d2 - d6)));
  }
  const double va = d3 * d6 - d5 * d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
    return distSq(p, addScaled(b, sub3(c, b),
                               (d4 - d3) / ((d4 - d3) + (d5 - d6))));
  }
  const double denom = 1.0 / (va + vb + vc);
  return distSq(p, addScaled(addScaled(a, ab, vb * denom), ac, vc * denom));
}

// Moller-Trumbore for a unit direction. A ray within 1e-12 radians of the
// triangle's plane is a miss here: if it met the triangle it also passed
// next to an edge, which the caller tests first.
bool rayHitsTriangle(const Vec3 &orig, const Vec3 &dir, const Vec3 &a,
                     const Vec3 &b, const Vec3 &c) {
  const auto e1 = sub3(b, a);
  const auto e2 = sub3(c, a);
  const auto pvec = crossP(dir, e2);
  const double det = dotP(e1, pvec);
  const auto n = crossP(e1, e2);
  if (std::fabs(det) <= 1e-12 * std::sqrt(dotP(n, n))) {
    return false;
  }
  const double inv = 1.0 / det;
  const auto tvec = sub3(orig, a);
  const double u = inv * dotP(tvec, pvec);
  if (u < 0.0 || u > 1.0) {
    return false;
  }
  const auto qvec = crossP(tvec, e1);
  const double v = inv * dotP(dir, qvec);
  if (v < 0.0 || u + v > 1.0) {
    return false;
  }
  return inv * dotP(e2, qvec) > 0.0;
}

struct Triangle {
  Vec3 a, b, c;
};

} // namespace

const std::vector<std::array<double, 3>> &insideRayDirections() {
  static const std::vector<std::array<double, 3>> dirs = [] {
    std::vector<std::array<double, 3>> raw = {
        {0.64278761, 0.58778525, 0.49240388},
        {-0.52991926, 0.71933980, 0.44964913},
        {0.37460659, -0.79863551, 0.47101007},
        {0.47942554, 0.33873792, -0.80934170},
        {-0.65605903, -0.42261826, 0.62524649},
        {0.84147098, -0.30901699, -0.44276253},
        {-0.29237170, 0.53582679, -0.79209936},
        {-0.41614684, -0.75470958, -0.50731958}};
    for (auto &d : raw) {
      const double inv = 1.0 / std::sqrt(dotP(d, d));
      d = {d[0] * inv, d[1] * inv, d[2] * inv};
    }
    return raw;
  }();
  return dirs;
}

InsideTrace
pointInFacesTrace(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                  const std::vector<std::vector<int>> &rings,
                  const std::vector<int> &faces, const std::array<double, 3> &guest,
                  const std::array<double, 3> *centreOpt) {
  InsideTrace out;
  std::vector<int> verts;
  for (int f : faces) {
    if (f < 0 || static_cast<std::size_t>(f) >= rings.size()) {
      continue;
    }
    verts.insert(verts.end(), rings[static_cast<std::size_t>(f)].begin(),
                 rings[static_cast<std::size_t>(f)].end());
  }
  std::sort(verts.begin(), verts.end());
  verts.erase(std::unique(verts.begin(), verts.end()), verts.end());
  if (verts.empty()) {
    return out;
  }
  const Vec3 centre = centreOpt != nullptr ? *centreOpt : periodicCentroid(yCloud, verts);
  auto unwrap = [&](double x, double y, double z) {
    const auto dr = minImage(yCloud, x, y, z, centre[0], centre[1], centre[2]);
    return Vec3{centre[0] + dr[0], centre[1] + dr[1], centre[2] + dr[2]};
  };

  // Fan each face from its own centroid, in the frame unwrapped about c_F.
  std::vector<Triangle> tris;
  double radius = 0.0;
  std::vector<Vec3> poly;
  for (int f : faces) {
    if (f < 0 || static_cast<std::size_t>(f) >= rings.size()) {
      continue;
    }
    const auto &ring = rings[static_cast<std::size_t>(f)];
    if (ring.size() < 3) {
      continue;
    }
    poly.clear();
    Vec3 m = {0.0, 0.0, 0.0};
    for (int atom : ring) {
      const auto &p = yCloud.pts[static_cast<std::size_t>(atom)];
      poly.push_back(unwrap(p.x, p.y, p.z));
      for (int k = 0; k < 3; k++) {
        m[static_cast<std::size_t>(k)] += poly.back()[static_cast<std::size_t>(k)];
      }
      radius = std::max(radius, std::sqrt(distSq(poly.back(), centre)));
    }
    const double inv = 1.0 / static_cast<double>(poly.size());
    m = {m[0] * inv, m[1] * inv, m[2] * inv};
    for (std::size_t k = 0; k < poly.size(); k++) {
      tris.push_back({m, poly[k], poly[(k + 1) % poly.size()]});
    }
  }
  if (tris.empty()) {
    return out;
  }

  const Vec3 g = unwrap(guest[0], guest[1], guest[2]);
  const double tol = kInsideRelTol * radius;
  const double tolSq = tol * tol;
  const double gc = std::sqrt(distSq(g, centre));
  // Every triangle lies in the ball of radius `radius` about c_F.
  if (gc > radius + tol) {
    return out;
  }
  for (const auto &t : tris) {
    if (triangleDistSq(g, t.a, t.b, t.c) <= tolSq) {
      out.onFace = true;
      out.inside = true;
      return out;
    }
  }

  // The segment g + s * dir, s in [0, reach], leaves that ball.
  const double reach = 2.0 * (gc + radius);
  const auto &dirs = insideRayDirections();
  auto clearanceSq = [&](const Vec3 &end, bool stopAtGraze) {
    double best = std::numeric_limits<double>::infinity();
    for (const auto &t : tris) {
      best = std::min({best, segmentDistSq(g, end, t.a, t.b),
                       segmentDistSq(g, end, t.b, t.c),
                       segmentDistSq(g, end, t.c, t.a)});
      if (stopAtGraze && best <= tolSq) {
        break;
      }
    }
    return best;
  };
  auto crossings = [&](const Vec3 &dir) {
    int hits = 0;
    for (const auto &t : tris) {
      hits += rayHitsTriangle(g, dir, t.a, t.b, t.c) ? 1 : 0;
    }
    return hits;
  };
  for (const auto &dir : dirs) {
    ++out.rays;
    if (clearanceSq(addScaled(g, dir, reach), true) <= tolSq) {
      continue;
    }
    out.crossings = crossings(dir);
    out.inside = (out.crossings % 2) == 1;
    return out;
  }
  // Every direction of the list grazes an edge: count along the one with
  // the widest clearance.
  std::size_t bestDir = 0;
  double bestClear = -1.0;
  for (std::size_t d = 0; d < dirs.size(); d++) {
    const double c = clearanceSq(addScaled(g, dirs[d], reach), false);
    if (c > bestClear) {
      bestClear = c;
      bestDir = d;
    }
  }
  out.crossings = crossings(dirs[bestDir]);
  out.inside = (out.crossings % 2) == 1;
  return out;
}

bool pointInFaces(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                  const std::vector<std::vector<int>> &rings,
                  const std::vector<int> &faces, const std::array<double, 3> &guest,
                  const std::array<double, 3> *centre) {
  return pointInFacesTrace(yCloud, rings, faces, guest, centre).inside;
}

GuestOccupancy
guestOccupancyInside(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                     const std::vector<std::vector<int>> &rings,
                     const std::vector<std::vector<int>> &cageFaces,
                     const std::vector<int> &guestIndices) {
  GuestOccupancy out;
  out.guestsPerCage.assign(cageFaces.size(), 0);
  std::vector<std::array<double, 3>> centres;
  centres.reserve(cageFaces.size());
  std::vector<std::vector<int>> cageVerts(cageFaces.size());
  for (std::size_t c = 0; c < cageFaces.size(); c++) {
    std::unordered_set<int> uniq;
    for (int f : cageFaces[c]) {
      if (f < 0 || static_cast<std::size_t>(f) >= rings.size()) {
        continue;
      }
      for (int v : rings[static_cast<std::size_t>(f)]) {
        uniq.insert(v);
      }
    }
    cageVerts[c].assign(uniq.begin(), uniq.end());
    std::sort(cageVerts[c].begin(), cageVerts[c].end());
    centres.push_back(periodicCentroid(yCloud, cageVerts[c]));
  }
  for (int g : guestIndices) {
    if (g < 0 || g >= yCloud.nop) {
      continue;
    }
    const auto &p = yCloud.pts[static_cast<std::size_t>(g)];
    const std::array<double, 3> gp = {p.x, p.y, p.z};
    int best = -1;
    double bestSq = std::numeric_limits<double>::infinity();
    for (std::size_t c = 0; c < cageFaces.size(); c++) {
      if (!pointInFaces(yCloud, rings, cageFaces[c], gp, &centres[c])) {
        continue;
      }
      const auto dr = minImage(yCloud, p.x, p.y, p.z, centres[c][0],
                               centres[c][1], centres[c][2]);
      const double d2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
      if (d2 <= bestSq) {
        bestSq = d2;
        best = static_cast<int>(c);
      }
    }
    out.cageOfGuest.push_back(best);
    out.centreDistance.push_back(best < 0 ? -1.0 : std::sqrt(bestSq));
    if (best < 0) {
      ++out.free;
    } else {
      ++out.guestsPerCage[static_cast<std::size_t>(best)];
    }
  }
  for (int n : out.guestsPerCage) {
    if (n > 0) {
      ++out.occupied;
    }
    if (n > 1) {
      ++out.multiply;
    }
  }
  int maxG = 0;
  for (int n : out.guestsPerCage) {
    maxG = std::max(maxG, n);
  }
  out.occupancyHistogram.assign(static_cast<std::size_t>(maxG) + 1, 0);
  for (int n : out.guestsPerCage) {
    ++out.occupancyHistogram[static_cast<std::size_t>(n)];
  }
  return out;
}

IceClusterIons
iceClusterIonCensus(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                    const std::vector<bool> &ice,
                    const std::vector<std::vector<int>> &nListByIndex,
                    const std::vector<int> &ionIndices, double cutoff) {
  IceClusterIons out;
  const int n = yCloud.nop;
  out.clusterOf.assign(static_cast<std::size_t>(n), -1);
  out.clusterOfIon.assign(ionIndices.size(), -1);
  std::vector<int> stack;
  for (int a = 0; a < n; a++) {
    if (a >= static_cast<int>(ice.size()) || !ice[static_cast<std::size_t>(a)] ||
        out.clusterOf[static_cast<std::size_t>(a)] >= 0) {
      continue;
    }
    const int cid = out.nClusters++;
    stack.push_back(a);
    out.clusterOf[static_cast<std::size_t>(a)] = cid;
    while (!stack.empty()) {
      const int u = stack.back();
      stack.pop_back();
      if (u < 0 || static_cast<std::size_t>(u) >= nListByIndex.size()) {
        continue;
      }
      for (std::size_t m = 1; m < nListByIndex[static_cast<std::size_t>(u)].size();
           m++) {
        const int v = nListByIndex[static_cast<std::size_t>(u)][m];
        if (v < 0 || v >= n || v >= static_cast<int>(ice.size())) {
          continue;
        }
        if (!ice[static_cast<std::size_t>(v)] ||
            out.clusterOf[static_cast<std::size_t>(v)] >= 0) {
          continue;
        }
        out.clusterOf[static_cast<std::size_t>(v)] = cid;
        stack.push_back(v);
      }
    }
  }
  out.ionsInCluster.assign(static_cast<std::size_t>(out.nClusters), 0);
  const double r2max = cutoff * cutoff;
  for (std::size_t k = 0; k < ionIndices.size(); k++) {
    const int g = ionIndices[k];
    if (g < 0 || g >= n) {
      continue;
    }
    int best = -1;
    double bestSq = r2max;
    for (int a = 0; a < n; a++) {
      if (out.clusterOf[static_cast<std::size_t>(a)] < 0) {
        continue;
      }
      const auto dr = gen::relDist(yCloud, g, a);
      const double d2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
      if (d2 <= bestSq) {
        bestSq = d2;
        best = a;
      }
    }
    if (best >= 0) {
      const int cid = out.clusterOf[static_cast<std::size_t>(best)];
      out.clusterOfIon[k] = cid;
      ++out.ionsInCluster[static_cast<std::size_t>(cid)];
    }
  }
  return out;
}

} // namespace site
