#include <water_select.hpp>

#include <generic.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace water {

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

bool isIonType(int type, const std::vector<int> &ionTypes) {
  return std::find(ionTypes.begin(), ionTypes.end(), type) != ionTypes.end();
}

bool typePresent(const Cloud &cloud, int type) {
  for (const auto &point : cloud.pts) {
    if (point.type == type) {
      return true;
    }
  }
  return false;
}

bool molPresent(const Cloud &cloud, int mol) {
  for (const auto &point : cloud.pts) {
    if (point.molID == mol) {
      return true;
    }
  }
  return false;
}

std::string missingType(int type) {
  return "atom type " + std::to_string(type) + " is not in the frame";
}

std::string missingMol(int mol) {
  return "molecule " + std::to_string(mol) + " is not in the frame";
}

std::vector<int> inferOxygenTypes(const Cloud &cloud,
                                  const std::vector<int> &ionTypes,
                                  const std::vector<int> &onlyMols) {
  const bool limit = !onlyMols.empty();
  std::unordered_set<int> allowed(onlyMols.begin(), onlyMols.end());
  std::unordered_map<int, std::vector<int>> groups;
  for (int i = 0; i < cloud.nop; ++i) {
    const int mol = cloud.pts[static_cast<std::size_t>(i)].molID;
    if (limit && allowed.find(mol) == allowed.end()) {
      continue;
    }
    groups[mol].push_back(i);
  }

  std::map<int, int> votes;
  int multi = 0;
  int ambiguous = 0;
  for (const auto &item : groups) {
    const auto &idx = item.second;
    if (idx.size() < 2) {
      continue;
    }
    ++multi;
    std::map<int, int> count;
    for (int i : idx) {
      const int type = cloud.pts[static_cast<std::size_t>(i)].type;
      if (isIonType(type, ionTypes)) {
        continue;
      }
      count[type] += 1;
    }
    std::vector<int> singles;
    for (const auto &entry : count) {
      if (entry.second == 1) {
        singles.push_back(entry.first);
      }
    }
    if (singles.size() == 1) {
      votes[singles[0]] += 1;
    } else if (singles.size() > 1) {
      ++ambiguous;
    }
  }
  if (multi > 0) {
    if (votes.size() == 1 && ambiguous == 0) {
      return {votes.begin()->first};
    }
    if (votes.size() > 1 || ambiguous > 0) {
      throw std::invalid_argument(
          "molecule has several singleton sites; pass --type");
    }
    throw std::invalid_argument("molecule has no water oxygen; pass --type");
  }

  std::map<int, int> typeCount;
  for (const auto &item : groups) {
    const int type =
        cloud.pts[static_cast<std::size_t>(item.second.front())].type;
    if (isIonType(type, ionTypes)) {
      continue;
    }
    typeCount[type] += static_cast<int>(item.second.size());
  }
  if (typeCount.size() == 1) {
    return {typeCount.begin()->first};
  }
  if (typeCount.empty()) {
    throw std::invalid_argument("selection keeps no water oxygens");
  }
  throw std::invalid_argument("mixed species need a water type or molecule");
}

double micDistance(const Cloud &cloud, double x, double y, double z, int j) {
  const auto &point = cloud.pts[static_cast<std::size_t>(j)];
  if (cloud.box.size() >= 6) {
    const auto dr =
        gen::triclinicMinImage(cloud, x, y, z, point.x, point.y, point.z);
    return std::sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
  }
  double delta[3] = {std::fabs(x - point.x), std::fabs(y - point.y),
                     std::fabs(z - point.z)};
  double r2 = 0.0;
  for (int k = 0; k < 3; ++k) {
    const double length = cloud.box[static_cast<std::size_t>(k)];
    if (length != 0.0) {
      delta[k] -= length * std::round(delta[k] / length);
    }
    r2 += delta[k] * delta[k];
  }
  return std::sqrt(r2);
}

} // namespace

bool iceLabel(molSys::atom_state_type label) {
  switch (label) {
  case molSys::atom_state_type::cubic:
  case molSys::atom_state_type::hexagonal:
  case molSys::atom_state_type::interfacial:
  case molSys::atom_state_type::clathrate:
  case molSys::atom_state_type::interClathrate:
  case molSys::atom_state_type::reCubic:
  case molSys::atom_state_type::reHex:
    return true;
  default:
    return false;
  }
}

bool hasMultiAtomMolecule(const Cloud &cloud) {
  std::unordered_map<int, int> count;
  for (int i = 0; i < cloud.nop; ++i) {
    const int n = ++count[cloud.pts[static_cast<std::size_t>(i)].molID];
    if (n >= 2) {
      return true;
    }
  }
  return false;
}

Selection selectOxygens(const Cloud &cloud, const Request &request) {
  std::vector<int> types = request.types;
  if (request.type > 0 &&
      std::find(types.begin(), types.end(), request.type) == types.end()) {
    types.insert(types.begin(), request.type);
  }
  for (int type : types) {
    if (!typePresent(cloud, type)) {
      throw std::invalid_argument(missingType(type));
    }
  }
  for (int type : types) {
    if (isIonType(type, request.ion_types)) {
      throw std::invalid_argument("atom type " + std::to_string(type) +
                                  " is an ion, not a water oxygen");
    }
  }
  for (int mol : request.molecules) {
    if (!molPresent(cloud, mol)) {
      throw std::invalid_argument(missingMol(mol));
    }
  }
  if (types.empty()) {
    if (request.mol_mode == MolMode::None) {
      throw std::invalid_argument("selection needs a water type or molecule");
    }
    const std::vector<int> limit = request.mol_mode == MolMode::Ids
                                       ? request.molecules
                                       : std::vector<int>{};
    types = inferOxygenTypes(cloud, request.ion_types, limit);
    for (int type : types) {
      if (!typePresent(cloud, type)) {
        throw std::invalid_argument(missingType(type));
      }
    }
  }

  const std::unordered_set<int> typeSet(types.begin(), types.end());
  const std::unordered_set<int> molSet(request.molecules.begin(),
                                       request.molecules.end());
  const bool molFilter = request.mol_mode == MolMode::Ids;
  Selection selection;
  selection.oxygen_types = types;
  selection.keep.assign(static_cast<std::size_t>(std::max(cloud.nop, 0)), 0);
  std::map<int, int> keptPerMol;
  for (int i = 0; i < cloud.nop; ++i) {
    const auto &point = cloud.pts[static_cast<std::size_t>(i)];
    if (typeSet.find(point.type) == typeSet.end()) {
      continue;
    }
    if (molFilter && molSet.find(point.molID) == molSet.end()) {
      continue;
    }
    selection.keep[static_cast<std::size_t>(i)] = 1;
    selection.index.push_back(i);
    keptPerMol[point.molID] += 1;
  }
  if (molFilter) {
    for (int mol : request.molecules) {
      if (keptPerMol[mol] == 0) {
        throw std::invalid_argument("molecule " + std::to_string(mol) +
                                    " has no water oxygen");
      }
    }
  }
  selection.n_oxygen = static_cast<int>(selection.index.size());
  if (selection.n_oxygen == 0) {
    throw std::invalid_argument("selection keeps no water oxygens");
  }
  return selection;
}

Cloud subsetCloud(const Cloud &cloud, const Selection &selection) {
  Cloud out;
  out.box = cloud.box;
  out.boxLow = cloud.boxLow;
  out.currentFrame = cloud.currentFrame;
  out.nop = selection.n_oxygen;
  out.pts.reserve(selection.index.size());
  for (int index : selection.index) {
    out.pts.push_back(cloud.pts[static_cast<std::size_t>(index)]);
    out.idIndexMap[out.pts.back().atomID] =
        static_cast<int>(out.pts.size()) - 1;
  }
  return out;
}

Cloud keepHistoricalOxygen(const Cloud &cloud) {
  bool has2 = false;
  bool has1 = false;
  for (const auto &point : cloud.pts) {
    has2 = has2 || point.type == 2;
    has1 = has1 || point.type == 1;
  }
  const int type = has2 ? 2 : 1;
  if (!has2 && !has1) {
    Cloud empty;
    empty.box = cloud.box;
    empty.boxLow = cloud.boxLow;
    empty.currentFrame = cloud.currentFrame;
    return empty;
  }
  Request request;
  request.type = type;
  return subsetCloud(cloud, selectOxygens(cloud, request));
}

std::vector<IonSite> collectIons(const Cloud &cloud,
                                 const Selection &selection) {
  std::unordered_set<int> waterMols;
  for (int index : selection.index) {
    waterMols.insert(cloud.pts[static_cast<std::size_t>(index)].molID);
  }
  std::vector<IonSite> ions;
  for (int i = 0; i < cloud.nop; ++i) {
    if (i < static_cast<int>(selection.keep.size()) &&
        selection.keep[static_cast<std::size_t>(i)] != 0) {
      continue;
    }
    const auto &point = cloud.pts[static_cast<std::size_t>(i)];
    if (waterMols.find(point.molID) != waterMols.end()) {
      continue;
    }
    ions.push_back(IonSite{point.x, point.y, point.z, point.type});
  }
  return ions;
}

IonReport ionFront(const Cloud &oxygen, const std::vector<IonSite> &ions) {
  IonReport report;
  if (ions.empty() || oxygen.nop <= 0 || oxygen.box.size() < 3) {
    return report;
  }
  report.any = true;
  double front = -1.0;
  bool sawIce = false;
  for (const IonSite &ion : ions) {
    int nearest = -1;
    double nearestDist = 0.0;
    for (int i = 0; i < oxygen.nop; ++i) {
      const double dist = micDistance(oxygen, ion.x, ion.y, ion.z, i);
      if (nearest < 0 || dist < nearestDist) {
        nearest = i;
        nearestDist = dist;
      }
      if (iceLabel(oxygen.pts[static_cast<std::size_t>(i)].iceType)) {
        if (!sawIce || dist < front) {
          front = dist;
          sawIce = true;
        }
      }
    }
    if (nearest < 0) {
      continue;
    }
    const auto label = oxygen.pts[static_cast<std::size_t>(nearest)].iceType;
    if (iceLabel(label)) {
      report.ice += 1;
    } else if (label == molSys::atom_state_type::water) {
      report.liquid += 1;
    }
  }
  if (sawIce) {
    report.front = front;
  }
  return report;
}

} // namespace water
