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

#include <chill_plus_scheme.hpp>

#include <bop.hpp>

#include <cctype>
#include <fstream>
#include <sstream>

namespace chill {
namespace {

std::string trim(std::string_view in) {
  std::size_t a = 0;
  std::size_t b = in.size();
  while (a < b && std::isspace(static_cast<unsigned char>(in[a]))) {
    ++a;
  }
  while (b > a && std::isspace(static_cast<unsigned char>(in[b - 1]))) {
    --b;
  }
  return std::string(in.substr(a, b - a));
}

std::string lower(std::string text) {
  for (char &c : text) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }
  return text;
}

std::string unquote(std::string text) {
  text = trim(text);
  if (text.size() >= 2) {
    const char q = text.front();
    if ((q == '"' || q == '\'') && text.back() == q) {
      text = text.substr(1, text.size() - 2);
    }
  }
  return trim(text);
}

ChillPlusRule namedRule(std::string name) {
  name = lower(unquote(name));
  if (name == "chill+" || name == "chill-plus" || name == "published") {
    return ChillPlusRule::Published;
  }
  if (name == "ovito-chill+" || name == "ovito-chill-plus") {
    return ChillPlusRule::Ovito;
  }
  throw std::invalid_argument(
      "unknown CHILL+ rule '" + name +
      "' (use chill+ or ovito-chill+)");
}

bool isKeyChar(unsigned char c) {
  return std::isalnum(c) || c == '_' || c == '-' || c == '+';
}

// The first `key` that is a YAML or Lua key. A `#` line comment, a `--`
// line comment, and a `--[[ ... ]]` long comment are not keys, and the
// word inside a quote is not a key. `rulex` is not `rule`.
std::string valueOfKey(const std::string &text, const std::string &key) {
  const std::string needle = lower(key);
  const std::size_t n = text.size();
  std::size_t i = 0;
  while (i < n) {
    if (text[i] == '#') {
      while (i < n && text[i] != '\n') {
        ++i;
      }
      continue;
    }
    if (text[i] == '-' && i + 1 < n && text[i + 1] == '-') {
      if (i + 3 < n && text[i + 2] == '[' && text[i + 3] == '[') {
        i += 4;
        while (i + 1 < n && !(text[i] == ']' && text[i + 1] == ']')) {
          ++i;
        }
        if (i + 1 < n) {
          i += 2;
        }
        continue;
      }
      while (i < n && text[i] != '\n') {
        ++i;
      }
      continue;
    }
    if (text[i] == '"' || text[i] == '\'') {
      const char quote = text[i];
      ++i;
      while (i < n && text[i] != quote && text[i] != '\n') {
        ++i;
      }
      if (i < n && text[i] == quote) {
        ++i;
      }
      continue;
    }
    if (isKeyChar(static_cast<unsigned char>(text[i])) &&
        (i == 0 || !isKeyChar(static_cast<unsigned char>(text[i - 1])))) {
      std::size_t end = i + 1;
      while (end < n && isKeyChar(static_cast<unsigned char>(text[end]))) {
        ++end;
      }
      if (lower(text.substr(i, end - i)) == needle) {
        std::size_t scan = end;
        while (scan < n &&
               std::isspace(static_cast<unsigned char>(text[scan]))) {
          ++scan;
        }
        if (scan < n && (text[scan] == ':' || text[scan] == '=')) {
          ++scan;
          while (scan < n &&
                 std::isspace(static_cast<unsigned char>(text[scan]))) {
            ++scan;
          }
          std::size_t stop = scan;
          while (stop < n && text[stop] != '\n' && text[stop] != ',' &&
                 text[stop] != '}' && text[stop] != '#') {
            ++stop;
          }
          return unquote(text.substr(scan, stop - scan));
        }
      }
      i = end;
      continue;
    }
    ++i;
  }
  return {};
}

} // namespace

int cutoffShellCount(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList, int iatom) {
  if (iatom < 0 || static_cast<std::size_t>(iatom) >= nList.size()) {
    return 0;
  }
  int shell = 0;
  const auto &row = nList[static_cast<std::size_t>(iatom)];
  for (std::size_t j = 1; j < row.size(); ++j) {
    if (yCloud.idIndexMap.find(row[j]) != yCloud.idIndexMap.end()) {
      ++shell;
    }
  }
  return shell;
}

ChillPlusRule chillPlusRule(std::string_view spec) {
  const std::string text = trim(spec);
  if (text.empty()) {
    return ChillPlusRule::Published;
  }
  const bool document = text.find('\n') != std::string::npos ||
                        text.find(':') != std::string::npos ||
                        text.find('=') != std::string::npos ||
                        text.find('{') != std::string::npos;
  if (!document) {
    return namedRule(text);
  }
  const std::string rule = valueOfKey(text, "rule");
  if (!rule.empty()) {
    return namedRule(rule);
  }
  const std::string named = valueOfKey(text, "chill-plus");
  if (!named.empty()) {
    return namedRule(named);
  }
  throw std::invalid_argument(
      "CHILL+ config has no rule (want rule: chill+ or rule: ovito-chill+)");
}

ChillPlusRule chillPlusRuleFromFile(const std::string &spec) {
  const std::string text = trim(spec);
  if (text.empty()) {
    return ChillPlusRule::Published;
  }
  std::ifstream in(text);
  if (!in) {
    return chillPlusRule(text);
  }
  std::ostringstream body;
  body << in.rdbuf();
  if (body.str().empty()) {
    throw std::invalid_argument("CHILL+ rule file is empty");
  }
  return chillPlusRule(body.str());
}

const char *chillPlusRuleName(ChillPlusRule rule) {
  return rule == ChillPlusRule::Ovito ? "ovito-chill+" : "chill+";
}

void applyChillPlus(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList, ChillPlusRule rule,
    bool isSlice) {
  if (rule == ChillPlusRule::Ovito) {
    // Whole cutoff row: q3 and c_ij both see the shell. Crystal rows
    // still require that row to hold four neighbours. OVITO's c_ij
    // edges are open.
    getCorrelPlus(yCloud, nList, isSlice, 0, false);
    getIceTypePlusNoPrint(yCloud, nList, isSlice, false);
    return;
  }
  getCorrelPlus(yCloud, nList, isSlice, 4);
  getIceTypePlusNoPrint(yCloud, nList, isSlice, true);
}

} // namespace chill
