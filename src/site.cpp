//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <site.hpp>

#include <generic.hpp>
#include <mol_sys.hpp>

#include <algorithm>
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
  IonEnvironment out;
  const int n = yCloud.nop;
  std::vector<char> isIon(static_cast<std::size_t>(std::max(n, 0)), 0);
  for (int i : ionIndices) {
    if (i >= 0 && i < n) {
      isIon[static_cast<std::size_t>(i)] = 1;
    }
  }
  const double cut2 = cutoff * cutoff;
  for (int i : ionIndices) {
    if (i < 0 || i >= n) {
      continue;
    }
    int shell = 0;
    int labelled = 0;
    for (int j = 0; j < n; j++) {
      if (j == i || isIon[static_cast<std::size_t>(j)]) {
        continue;
      }
      if (waterType != 0 && yCloud.pts[static_cast<std::size_t>(j)].type != waterType) {
        continue;
      }
      if (gen::periodicDistSq(yCloud, i, j) >= cut2) {
        continue;
      }
      ++shell;
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
  }
  return out;
}

} // namespace site
