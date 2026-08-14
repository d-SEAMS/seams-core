//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#ifndef SEAMS_YAML_H_
#define SEAMS_YAML_H_

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

/** @file seams_yaml.hpp
 *  @brief Tiny reader for the yodaStruct config schema.
 */

namespace seams {

struct LuaConfig {
  std::string trajectory;
  std::string variables;
  bool bulkUse = false;
  bool bulkTopo = false;
  bool bulkBop = false;
  bool topoOneDim = false;
  bool topoTwoDim = false;
  bool structureDesc = false;
};

inline std::string trim(const std::string &s) {
  const auto a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) {
    return {};
  }
  const auto b = s.find_last_not_of(" \t\r\n");
  return s.substr(a, b - a + 1);
}

inline std::string unquote(std::string s) {
  if (s.size() >= 2 && ((s.front() == '"' && s.back() == '"') ||
                        (s.front() == '\'' && s.back() == '\''))) {
    return s.substr(1, s.size() - 2);
  }
  return s;
}

inline bool asBool(const std::string &v) {
  return v == "true" || v == "True" || v == "yes" || v == "1";
}

inline LuaConfig loadLuaConfig(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("cannot open config " + path);
  }
  LuaConfig cfg;
  std::string section;
  std::string line;
  while (std::getline(in, line)) {
    const auto hash = line.find('#');
    if (hash != std::string::npos) {
      line = line.substr(0, hash);
    }
    line = trim(line);
    if (line.empty()) {
      continue;
    }
    if (line.back() == ':') {
      section = trim(line.substr(0, line.size() - 1));
      continue;
    }
    const auto colon = line.find(':');
    if (colon == std::string::npos) {
      continue;
    }
    const std::string key = trim(line.substr(0, colon));
    const std::string val = unquote(trim(line.substr(colon + 1)));
    if (section.empty()) {
      if (key == "trajectory") {
        cfg.trajectory = val;
      } else if (key == "variables") {
        cfg.variables = val;
      }
      continue;
    }
    if (section == "bulk") {
      if (key == "use") {
        cfg.bulkUse = asBool(val);
      } else if (key == "topologicalNetworkCriterion") {
        cfg.bulkTopo = asBool(val);
      } else if (key == "bondOrderParameters") {
        cfg.bulkBop = asBool(val);
      }
    } else if (section == "topoOneDim" && key == "use") {
      cfg.topoOneDim = asBool(val);
    } else if (section == "topoTwoDim" && key == "use") {
      cfg.topoTwoDim = asBool(val);
    } else if (section == "structureDesc" && key == "use") {
      cfg.structureDesc = asBool(val);
    }
  }
  return cfg;
}

} // namespace seams

#endif // SEAMS_YAML_H_
