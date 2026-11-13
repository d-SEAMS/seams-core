//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//-----------------------------------------------------------------------------------

#include <bop.hpp>
#include <chill_offload.hpp>
#include <chill_plus_scheme.hpp>
#include <ozaki.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_c_api.h>
#include <seams_input.hpp>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <type_traits>
#include <utility>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

struct Frame {
  Cloud cloud;
  int focus = 0;
};

Frame readXyz(const std::string &path) {
  std::ifstream in(path);
  REQUIRE(in);
  Frame frame;
  int n = 0;
  in >> n;
  double lx = 0, ly = 0, lz = 0;
  in >> lx >> ly >> lz;
  frame.cloud.nop = n;
  frame.cloud.box = {lx, ly, lz};
  frame.cloud.boxLow = {0, 0, 0};
  frame.cloud.pts.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    double x = 0, y = 0, z = 0;
    in >> x >> y >> z;
    auto &p = frame.cloud.pts[static_cast<std::size_t>(i)];
    p.x = x;
    p.y = y;
    p.z = z;
    p.type = 1;
    p.atomID = i + 1;
    p.molID = i + 1;
    frame.cloud.idIndexMap[p.atomID] = i;
  }
  REQUIRE(static_cast<bool>(in >> frame.focus));
  return frame;
}

Cloud diamond(int nrep) {
  const double oo = 2.75;
  const double a = oo * 4.0 / std::sqrt(3.0);
  const double basis[2][3] = {{0, 0, 0}, {0.25, 0.25, 0.25}};
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  Cloud cloud;
  cloud.box = {nrep * a, nrep * a, nrep * a};
  cloud.boxLow = {0, 0, 0};
  for (int i = 0; i < nrep; ++i) {
    for (int j = 0; j < nrep; ++j) {
      for (int k = 0; k < nrep; ++k) {
        for (const auto &shift : fcc) {
          for (const auto &site : basis) {
            molSys::Point<double> p;
            p.x = (i + shift[0] + site[0]) * a;
            p.y = (j + shift[1] + site[1]) * a;
            p.z = (k + shift[2] + site[2]) * a;
            p.type = 1;
            p.atomID = cloud.nop + 1;
            p.molID = p.atomID;
            cloud.idIndexMap[p.atomID] = cloud.nop;
            cloud.pts.push_back(p);
            cloud.nop += 1;
          }
        }
      }
    }
  }
  return cloud;
}

int labelAt(const Cloud &cloud, int i) {
  return static_cast<int>(cloud.pts[static_cast<std::size_t>(i)].iceType);
}

void requireUniform(const std::string &path, int type, molSys::atom_state_type want) {
  Cloud cloud;
  cloud = sinp::readLammpsTrj(path, 1, cloud);
  REQUIRE(cloud.nop > 0);
  auto rows = nneigh::neighListO(3.5, cloud, type);
  auto ovito = cloud;
  chill::applyChillPlus(cloud, rows, chill::ChillPlusRule::Published, false);
  chill::applyChillPlus(ovito, rows, chill::ChillPlusRule::Ovito, false);
  for (int i = 0; i < cloud.nop; ++i) {
    if (cloud.pts[static_cast<std::size_t>(i)].type != type) {
      continue;
    }
    INFO(path << " atom " << i);
    REQUIRE(cloud.pts[static_cast<std::size_t>(i)].iceType == want);
    REQUIRE(ovito.pts[static_cast<std::size_t>(i)].iceType == want);
  }
}

} // namespace

TEST_CASE("CHILL+ rule text accepts the two names and rejects the rest",
          "[chill][rule]") {
  REQUIRE(chill::chillPlusRule("") == chill::ChillPlusRule::Published);
  REQUIRE(chill::chillPlusRule("  chill+ ") == chill::ChillPlusRule::Published);
  REQUIRE(chill::chillPlusRule("CHILL-PLUS") == chill::ChillPlusRule::Published);
  REQUIRE(chill::chillPlusRule("published") == chill::ChillPlusRule::Published);
  REQUIRE(chill::chillPlusRule("ovito-chill+") == chill::ChillPlusRule::Ovito);
  REQUIRE(chill::chillPlusRule("OVITO-CHILL-PLUS") == chill::ChillPlusRule::Ovito);
  REQUIRE(std::string(chill::chillPlusRuleName(chill::ChillPlusRule::Published)) ==
          "chill+");
  REQUIRE(std::string(chill::chillPlusRuleName(chill::ChillPlusRule::Ovito)) ==
          "ovito-chill+");
  REQUIRE_THROWS_AS(chill::chillPlusRule("chill"), std::invalid_argument);
  REQUIRE_THROWS_AS(chill::chillPlusRule("ovito"), std::invalid_argument);
  REQUIRE_THROWS_AS(chill::chillPlusRule("rule:"), std::invalid_argument);
  REQUIRE(chill::chillPlusRule("rule: ovito-chill+") == chill::ChillPlusRule::Ovito);
  REQUIRE(chill::chillPlusRule("chill-plus: chill+") == chill::ChillPlusRule::Published);
  REQUIRE(chill::chillPlusRule("{rule = \"ovito-chill+\"}") ==
          chill::ChillPlusRule::Ovito);
  REQUIRE(chill::chillPlusRule("rule='chill+'") == chill::ChillPlusRule::Published);
  REQUIRE_THROWS_AS(chill::chillPlusRule("{cutoff = 3.5}"), std::invalid_argument);
  REQUIRE_THROWS_AS(chill::chillPlusRule("foo: rule"), std::invalid_argument);
  REQUIRE_THROWS_AS(chill::chillPlusRule("rulex:"), std::invalid_argument);
  REQUIRE(chill::chillPlusRule("rule: ovito-chill+\n") ==
          chill::ChillPlusRule::Ovito);
  REQUIRE(chill::chillPlusRule("rule: ovito-chill+, extra") ==
          chill::ChillPlusRule::Ovito);
  REQUIRE_THROWS_AS(chill::chillPlusRule("{ovito}"), std::invalid_argument);
  REQUIRE(chill::chillPlusRule(
              "-- the rule option is not a key\nrule = \"ovito-chill+\"") ==
          chill::ChillPlusRule::Ovito);
  REQUIRE(chill::chillPlusRule(
              "-- rule = \"ovito-chill+\"\nrule = \"chill+\"") ==
          chill::ChillPlusRule::Published);
  REQUIRE(chill::chillPlusRule(
              "# rule: ovito-chill+\nrule: chill+") ==
          chill::ChillPlusRule::Published);
  REQUIRE(chill::chillPlusRule(
              "note = \"rule: ovito-chill+\"\nrule: chill+") ==
          chill::ChillPlusRule::Published);
  REQUIRE(chill::chillPlusRule(
              "--[[ rule = \"ovito-chill+\" ]]\nrule = \"chill+\"") ==
          chill::ChillPlusRule::Published);
  REQUIRE_THROWS_AS(
      chill::chillPlusRule("--[[ rule = \"ovito-chill+\""),
      std::invalid_argument);
  REQUIRE(chill::chillPlusRuleFromFile("../examples/chill-plus.lua") ==
          chill::ChillPlusRule::Published);
  REQUIRE(chill::chillPlusRuleFromFile("../examples/chill-plus.yaml") ==
          chill::ChillPlusRule::Published);
  REQUIRE(chill::chillPlusRuleFromFile("../examples/ovito-chill-plus.yaml") ==
          chill::ChillPlusRule::Ovito);
  REQUIRE(chill::chillPlusRuleFromFile("chill+") == chill::ChillPlusRule::Published);
  REQUIRE(chill::chillPlusRuleFromFile("ovito-chill+") == chill::ChillPlusRule::Ovito);
  REQUIRE(chill::chillPlusRule("rule: chill+#") == chill::ChillPlusRule::Published);
  REQUIRE(chill::chillPlusRule("--\nrule: chill+") == chill::ChillPlusRule::Published);
  REQUIRE(chill::chillPlusRule("--[[rule]]\nrule: chill+") ==
          chill::ChillPlusRule::Published);
  REQUIRE_THROWS_AS(chill::chillPlusRule("rule: ''"), std::invalid_argument);
  REQUIRE(chill::chillPlusRule("-\nrule = chill+") == chill::ChillPlusRule::Published);
  REQUIRE(chill::chillPlusRule("-\nrule = chill+") == chill::ChillPlusRule::Published);
  REQUIRE(chill::chillPlusRule("cutoff: -1\nrule: chill+") ==
          chill::ChillPlusRule::Published);
  REQUIRE_THROWS_AS(chill::chillPlusRule("rule: \"chill+"), std::invalid_argument);
  REQUIRE(chill::chillPlusRule("--[ not a long comment\nrule: chill+") ==
          chill::ChillPlusRule::Published);
  REQUIRE(chill::chillPlusRule("rule: chill+\n#") ==
          chill::ChillPlusRule::Published);
  REQUIRE_THROWS_AS(chill::chillPlusRule("-- note"), std::invalid_argument);
  REQUIRE(chill::chillPlusRule("#\nrule: chill+") == chill::ChillPlusRule::Published);
  REQUIRE(chill::chillPlusRule("-\nrule: chill+") == chill::ChillPlusRule::Published);
  REQUIRE(chill::chillPlusRule("--[\nrule: chill+") == chill::ChillPlusRule::Published);
  REQUIRE_THROWS_AS(chill::chillPlusRule("rule: \"\""), std::invalid_argument);
  REQUIRE(chill::chillPlusRule("rule = 'chill+'") == chill::ChillPlusRule::Published);
  REQUIRE_THROWS_AS(chill::chillPlusRule("note = \"unterminated"),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(chill::chillPlusRule("{\n# no newline"), std::invalid_argument);
  REQUIRE_THROWS_AS(chill::chillPlusRule("{\n-"), std::invalid_argument);
  REQUIRE_THROWS_AS(chill::chillPlusRule("{\n--"), std::invalid_argument);
  REQUIRE(chill::chillPlusRule("--[[ ]x ]]\nrule: chill+") ==
          chill::ChillPlusRule::Published);
  REQUIRE_THROWS_AS(chill::chillPlusRule("{\n-- note"), std::invalid_argument);
  REQUIRE(chill::chillPlusRule("note: \"hi\nrule: chill+") ==
          chill::ChillPlusRule::Published);
  REQUIRE(chill::chillPlusRule("\nrule: chill+") == chill::ChillPlusRule::Published);
  REQUIRE_THROWS_AS(chill::chillPlusRule("rule chill+\n"), std::invalid_argument);
  REQUIRE(chill::chillPlusRule("rule=\nchill+") == chill::ChillPlusRule::Published);
  REQUIRE(chill::cutoffShellCount(Cloud{}, {}, -1) == 0);
  REQUIRE(chill::cutoffShellCount(Cloud{}, {}, 3) == 0);
  Cloud ghost;
  ghost.nop = 1;
  ghost.pts.resize(1);
  ghost.pts[0].atomID = 1;
  ghost.idIndexMap[1] = 0;
  ghost.pts[0].type = 1;
  std::vector<std::vector<int>> rows{{1, 99}};
  REQUIRE(chill::cutoffShellCount(ghost, rows, 0) == 0);
}

TEST_CASE("a CHILL+ rule file is read, and an empty file fails", "[chill][rule]") {
  const auto dir = std::filesystem::temp_directory_path() / "dseams-chill-rule";
  std::filesystem::create_directories(dir);
  const auto path = dir / "rule-ovito.yaml";
  {
    std::ofstream out(path);
    out << "rule: ovito-chill+\n";
  }
  REQUIRE(chill::chillPlusRuleFromFile(path.string()) == chill::ChillPlusRule::Ovito);
  const auto empty = dir / "rule-empty.yaml";
  {
    std::ofstream out(empty);
  }
  REQUIRE_THROWS_AS(chill::chillPlusRuleFromFile(empty.string()),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(chill::chillPlusRuleFromFile("chill-def/no-such-rule"),
                    std::invalid_argument);
  REQUIRE(chill::chillPlusRuleFromFile("") == chill::ChillPlusRule::Published);
  std::filesystem::remove_all(dir);
}

TEST_CASE("the shell gate sends a five-neighbour oxygen to water", "[chill][rule]") {
  auto frame = readXyz("chill-def/counterexample-five.txt");
  auto rows = nneigh::neighListO(3.5, frame.cloud, 1);
  REQUIRE(chill::cutoffShellCount(frame.cloud, rows, frame.focus) == 5);
  auto published = frame.cloud;
  auto ovito = frame.cloud;
  chill::applyChillPlus(published, rows, chill::ChillPlusRule::Published, false);
  chill::applyChillPlus(ovito, rows, chill::ChillPlusRule::Ovito, false);
  REQUIRE(labelAt(published, frame.focus) ==
          static_cast<int>(molSys::atom_state_type::water));
  REQUIRE(labelAt(ovito, frame.focus) ==
          static_cast<int>(molSys::atom_state_type::water));
  int cubic = 0;
  for (int i = 0; i < published.nop; ++i) {
    if (published.pts[static_cast<std::size_t>(i)].iceType ==
        molSys::atom_state_type::cubic) {
      ++cubic;
    }
  }
  REQUIRE(cubic > 0);
}

TEST_CASE("q3 over the shell and the neighbour test are what ovito-chill+ changes",
          "[chill][rule]") {
  auto q3 = readXyz("chill-def/counterexample-q3.txt");
  auto rows = nneigh::neighListO(3.5, q3.cloud, 1);
  REQUIRE(chill::cutoffShellCount(q3.cloud, rows, q3.focus) == 4);
  auto published = q3.cloud;
  auto ovito = q3.cloud;
  chill::applyChillPlus(published, rows, chill::ChillPlusRule::Published, false);
  chill::applyChillPlus(ovito, rows, chill::ChillPlusRule::Ovito, false);
  REQUIRE(labelAt(published, q3.focus) ==
          static_cast<int>(molSys::atom_state_type::water));
  REQUIRE(labelAt(ovito, q3.focus) ==
          static_cast<int>(molSys::atom_state_type::interfacial));

  auto face = readXyz("chill-def/counterexample-interface.txt");
  rows = nneigh::neighListO(3.5, face.cloud, 1);
  REQUIRE(chill::cutoffShellCount(face.cloud, rows, face.focus) == 4);
  published = face.cloud;
  ovito = face.cloud;
  chill::applyChillPlus(published, rows, chill::ChillPlusRule::Published, false);
  chill::applyChillPlus(ovito, rows, chill::ChillPlusRule::Ovito, false);
  REQUIRE(labelAt(published, face.focus) ==
          static_cast<int>(molSys::atom_state_type::water));
  REQUIRE(labelAt(ovito, face.focus) ==
          static_cast<int>(molSys::atom_state_type::interfacial));
}

TEST_CASE("plain CHILL still scores four nearest bonds when a fifth oxygen is inside 3.5",
          "[chill][rule]") {
  auto cloud = diamond(2);
  molSys::Point<double> extra;
  extra.x = cloud.pts[0].x + 3.4;
  extra.y = cloud.pts[0].y;
  extra.z = cloud.pts[0].z;
  extra.type = 1;
  extra.atomID = cloud.nop + 1;
  extra.molID = extra.atomID;
  cloud.idIndexMap[extra.atomID] = cloud.nop;
  cloud.pts.push_back(extra);
  cloud.nop += 1;
  auto rows = nneigh::neighListO(3.5, cloud, 1);
  REQUIRE(chill::cutoffShellCount(cloud, rows, 0) == 5);
  chill::getCorrel(cloud, rows, false);
  chill::getIceTypeNoPrint(cloud, rows, false);
  REQUIRE(cloud.pts[0].iceType == molSys::atom_state_type::cubic);
  chill::applyChillPlus(cloud, rows, chill::ChillPlusRule::Published, false);
  REQUIRE(cloud.pts[0].iceType == molSys::atom_state_type::water);
}

TEST_CASE("ideal Ih, Ic, sI and sII are the same crystal under both CHILL+ rules",
          "[chill][rule][fixture]") {
  requireUniform("traj/ice1h.lammpstrj", 1, molSys::atom_state_type::hexagonal);
  requireUniform("traj/ice1h.lammpstrj", 2, molSys::atom_state_type::hexagonal);
  requireUniform("traj/genice_sI.lammpstrj", 1, molSys::atom_state_type::clathrate);
  requireUniform("traj/genice_sII.lammpstrj", 1, molSys::atom_state_type::clathrate);
  Cloud cubic = diamond(2);
  auto rows = nneigh::neighListO(3.5, cubic, 1);
  chill::applyChillPlus(cubic, rows, chill::ChillPlusRule::Published, false);
  for (const auto &point : cubic.pts) {
    REQUIRE(point.iceType == molSys::atom_state_type::cubic);
  }
}

TEST_CASE("the C API rejects an unknown CHILL+ rule and keeps the published default",
          "[chill][rule]") {
  auto cloud = diamond(2);
  std::vector<double> xyz(static_cast<std::size_t>(cloud.nop) * 3);
  for (int i = 0; i < cloud.nop; ++i) {
    xyz[static_cast<std::size_t>(3 * i)] = cloud.pts[static_cast<std::size_t>(i)].x;
    xyz[static_cast<std::size_t>(3 * i + 1)] = cloud.pts[static_cast<std::size_t>(i)].y;
    xyz[static_cast<std::size_t>(3 * i + 2)] = cloud.pts[static_cast<std::size_t>(i)].z;
  }
  const double box[6] = {0, 0, 0, cloud.box[0], cloud.box[1], cloud.box[2]};
  std::vector<int> labels(static_cast<std::size_t>(cloud.nop), -1);
  REQUIRE(seams_chill_plus_rule(xyz.data(), cloud.nop, box, 3.5, "cutoff",
                                "not-a-rule", labels.data()) == 1);
  REQUIRE(seams_chill_plus_flat(nullptr, 0, nullptr, nullptr) == 1);
  REQUIRE(seams_chill_plus_flat(xyz.data(), cloud.nop, box, labels.data()) == 0);
  for (int label : labels) {
    REQUIRE(label == static_cast<int>(molSys::atom_state_type::cubic));
  }
  REQUIRE(seams_chill_plus_rule(xyz.data(), cloud.nop, box, 3.5, "cutoff",
                                "ovito-chill+", labels.data()) == 0);
  for (int label : labels) {
    REQUIRE(label == static_cast<int>(molSys::atom_state_type::cubic));
  }
}

TEST_CASE("the flattened CHILL+ path uses the published shell count", "[chill][rule]") {
  Cloud iceCloud = diamond(2);
  auto rows = nneigh::neighListO(3.5, iceCloud, 1);
  const auto flat = chill::specializedChillPlus(iceCloud, rows);
  REQUIRE(flat.iceType.size() == static_cast<std::size_t>(iceCloud.nop));
  for (int label : flat.iceType) {
    REQUIRE(label == static_cast<int>(molSys::atom_state_type::cubic));
  }
  auto five = readXyz("chill-def/counterexample-five.txt");
  auto fiveRows = nneigh::neighListO(3.5, five.cloud, 1);
  fiveRows[0].push_back(999999);
  const auto fiveFlat = chill::specializedChillPlus(five.cloud, fiveRows);
  REQUIRE(fiveFlat.iceType[0] == static_cast<int>(molSys::atom_state_type::water));
  auto face = readXyz("chill-def/counterexample-interface.txt");
  auto faceRows = nneigh::neighListO(3.5, face.cloud, 1);
  const auto faceFlat = chill::specializedChillPlus(face.cloud, faceRows);
  REQUIRE(faceFlat.iceType.size() == static_cast<std::size_t>(face.cloud.nop));
  auto q3 = readXyz("chill-def/counterexample-q3.txt");
  auto q3Rows = nneigh::neighListO(3.5, q3.cloud, 1);
  const auto q3Flat = chill::specializedChillPlus(q3.cloud, q3Rows);
  REQUIRE(q3Flat.iceType.size() == static_cast<std::size_t>(q3.cloud.nop));
  Cloud empty;
  const auto none = chill::specializedChillPlus(empty, {});
  REQUIRE(none.iceType.empty());
  unsetenv("SEAMS_OFFLOAD");
  REQUIRE(chill::preferOffload() == false);
  const auto host = chill::chillPlus(iceCloud, rows);
  REQUIRE(host.usedDevice == false);
  REQUIRE(host.iceType[0] == static_cast<int>(molSys::atom_state_type::cubic));
  setenv("SEAMS_OFFLOAD", "1", 1);
  REQUIRE(chill::preferOffload() == true);
  const auto again = chill::chillPlus(iceCloud, rows);
  REQUIRE(again.iceType[0] == static_cast<int>(molSys::atom_state_type::cubic));
  setenv("SEAMS_OFFLOAD", "0", 1);
  REQUIRE(chill::preferOffload() == false);
  unsetenv("SEAMS_OFFLOAD");
  for (auto &pt : iceCloud.pts) {
    pt.inSlice = false;
  }
  iceCloud.pts[0].inSlice = true;
  chill::applyChillPlus(iceCloud, rows, chill::ChillPlusRule::Published, true);
  REQUIRE(iceCloud.pts[0].iceType == molSys::atom_state_type::cubic);
  chill::applyChillPlus(iceCloud, rows, chill::ChillPlusRule::Ovito, true);
  REQUIRE(iceCloud.pts[0].iceType == molSys::atom_state_type::cubic);
}

Cloud boundCloud() {
  Cloud cloud;
  cloud.box = {20, 20, 20};
  cloud.boxLow = {0, 0, 0};
  for (int i = 0; i < 5; ++i) {
    molSys::Point<double> point;
    point.x = i == 0 ? 0 : 1.5;
    point.y = i == 0 ? 0 : (i - 1);
    point.z = 0;
    point.type = 1;
    point.atomID = i + 1;
    point.molID = point.atomID;
    cloud.idIndexMap[point.atomID] = cloud.nop;
    cloud.pts.push_back(point);
    cloud.nop += 1;
  }
  return cloud;
}

void stampBonds(Cloud &cloud, double cij, bool closed) {
  chill::BondClassifier rule = chill::chillPlusRule();
  cloud.pts[0].c_ij.clear();
  for (int bond = 0; bond < 4; ++bond) {
    molSys::Result result;
    result.c_value = cij;
    result.classifier = chill::classifyBond(cij, rule, closed);
    cloud.pts[0].c_ij.push_back(result);
  }
}

TEST_CASE("a 2.13 BondClassifier passed by value keeps inclusive edges",
          "[chill][abi]") {
  chill::BondClassifier rule = chill::chillPlusRule();
  static_assert(std::is_standard_layout_v<chill::BondClassifier>);
  static_assert(offsetof(chill::BondClassifier, coordinationNumber) == 24);
  alignas(chill::BondClassifier) unsigned char raw[sizeof(rule)];
  std::memcpy(raw, &rule, sizeof(rule));
  // A 2.13 caller writes three doubles and one int. The tail of the
  // 32-byte value is padding. Zero it, which is what a by-value pass
  // leaves behind when that byte is not a member.
  REQUIRE(sizeof(rule) > 28);
  std::memset(raw + 28, 0, sizeof(rule) - 28);
  std::memcpy(&rule, raw, sizeof(rule));
  REQUIRE(chill::classifyBond(-0.8, rule) == molSys::bond_type::staggered);
  REQUIRE(chill::classifyBond(-0.35, rule) == molSys::bond_type::eclipsed);
  REQUIRE(chill::classifyBond(0.25, rule) == molSys::bond_type::eclipsed);
  REQUIRE(chill::classifyBond(-0.8, rule, false) == molSys::bond_type::out_of_range);
}

TEST_CASE("a correlation exactly on a bound follows the edge each rule uses",
          "[chill][rule]") {
  chill::BondClassifier rule = chill::chillPlusRule();
  REQUIRE(chill::classifyBond(-0.8, rule, true) == molSys::bond_type::staggered);
  REQUIRE(chill::classifyBond(-0.8, rule, false) == molSys::bond_type::out_of_range);
  REQUIRE(chill::classifyBond(-0.35, rule, true) == molSys::bond_type::eclipsed);
  REQUIRE(chill::classifyBond(-0.35, rule, false) == molSys::bond_type::out_of_range);
  REQUIRE(chill::classifyBond(0.25, rule, true) == molSys::bond_type::eclipsed);
  REQUIRE(chill::classifyBond(0.25, rule, false) == molSys::bond_type::out_of_range);
  REQUIRE(chill::classifyBond(-1.0, rule, false) == molSys::bond_type::staggered);
  REQUIRE(chill::classifyBond(0.0, rule, false) == molSys::bond_type::eclipsed);
  REQUIRE(chill::classifyBond(0.5, rule, true) == molSys::bond_type::out_of_range);

  const std::vector<std::vector<int>> rows{{1, 2, 3, 4, 5}};
  Cloud onStaggered = boundCloud();
  stampBonds(onStaggered, -0.8, true);
  chill::getIceTypePlusNoPrint(onStaggered, rows, false, true);
  REQUIRE(onStaggered.pts[0].iceType == molSys::atom_state_type::cubic);
  stampBonds(onStaggered, -0.8, false);
  chill::getIceTypePlusNoPrint(onStaggered, rows, false, false);
  REQUIRE(onStaggered.pts[0].iceType == molSys::atom_state_type::water);

  const auto dir = std::filesystem::temp_directory_path() / "chill-plus-bop";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directory(dir);
  Cloud printed = diamond(1);
  auto printedRows = nneigh::neighListO(3.5, printed, 1);
  chill::getCorrelPlus(printed, printedRows, false, 4, true);
  chill::getIceTypePlus(printed, printedRows, dir.string() + "/", 1, false,
                        "chillPlus.txt", true);
  REQUIRE(std::filesystem::exists(dir / "bop" / "chillPlus.txt"));
  const auto saved = seams::product::reduction();
  seams::product::setReduction(seams::product::Reduction::Pairwise);
  chill::getCorrelPlus(printed, printedRows, false, 4, false);
  seams::product::setReduction(saved);
  REQUIRE(printed.pts[0].c_ij.size() == 4);
}

TEST_CASE("CHILL, q6 and the C API error returns are exercised", "[chill][rule]") {
  Cloud cloud = diamond(2);
  auto rows = nneigh::neighListO(3.5, cloud, 1);
  chill::getCorrel(cloud, rows, false);
  const auto dir = std::filesystem::temp_directory_path() / "chill-plain";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directory(dir);
  chill::getIceType(cloud, rows, dir.string() + "/", 1, false, "chill.txt");
  REQUIRE(std::filesystem::exists(dir / "bop" / "chill.txt"));
  const auto q6 = chill::getq6(cloud, rows, false);
  REQUIRE(q6.size() == static_cast<std::size_t>(cloud.nop));
  REQUIRE_FALSE(chill::bondClassifierNames().empty());
  REQUIRE(chill::printIceType(cloud, dir.string() + "/", 1) >= 0);
  auto nudged = q6;
  chill::reclassifyWater(cloud, nudged);
  const std::array<double, 2> angles{{0.2, 0.4}};
  REQUIRE(sph::lookupTableQ3(-1, angles) == std::complex<double>(0.0, 0.0));
  REQUIRE(sph::lookupTableQ6(99, angles) == std::complex<double>(0.0, 0.0));
  REQUIRE(sph::lookupTableQ4(-1, angles) == std::complex<double>(0.0, 0.0));
  REQUIRE(sph::lookupTableQ8(-1, angles) == std::complex<double>(0.0, 0.0));

  Cloud ice;
  ice = sinp::readLammpsTrj("traj/ice1h.lammpstrj", 1, ice);
  auto iceRows = nneigh::neighListO(3.5, ice, 1);
  const auto flat = chill::specializedChillPlus(ice, iceRows);
  REQUIRE(flat.iceType.size() == static_cast<std::size_t>(ice.nop));
  setenv("SEAMS_OFFLOAD", "", 1);
  REQUIRE(chill::preferOffload() == false);
  setenv("SEAMS_OFFLOAD", "0", 1);
  REQUIRE(chill::preferOffload() == false);
  setenv("SEAMS_OFFLOAD", "00", 1);
  REQUIRE(chill::preferOffload() == true);
  setenv("SEAMS_OFFLOAD", "2", 1);
  REQUIRE(chill::preferOffload() == true);
  unsetenv("SEAMS_OFFLOAD");

  double xyz[3] = {0, 0, 0};
  double box[6] = {0, 0, 0, 10, 10, 10};
  int label = -1;
  REQUIRE(seams_chill_plus_graph(nullptr, 1, box, 3.5, "cutoff", &label) == 1);
  REQUIRE(seams_chill_plus_graph(xyz, 1, nullptr, 3.5, "cutoff", &label) == 1);
  REQUIRE(seams_chill_plus_graph(xyz, 1, box, 3.5, "cutoff", nullptr) == 1);
  REQUIRE(seams_chill_plus_graph(xyz, 0, box, 3.5, "cutoff", &label) == 1);
  REQUIRE(seams_chill_plus_graph(xyz, 1, box, 0.0, "cutoff", &label) == 1);
  REQUIRE(seams_chill_plus_graph(xyz, 1, box, 3.5, nullptr, &label) == 0);
  REQUIRE(seams_chill_plus_graph(xyz, 1, box, 3.5, "", &label) == 0);
  REQUIRE(seams_chill_plus_graph(xyz, 1, box, 3.5, "nope", &label) == 1);
  REQUIRE(seams_chill_plus_graph(xyz, 1, box, 3.5, "mutual", &label) == 0);
  double thin[6] = {0, 0, 0, 0, 10, 10};
  REQUIRE(seams_chill_plus_graph(xyz, 1, thin, 3.5, "cutoff", &label) == 1);
  thin[3] = 10;
  thin[4] = 0;
  REQUIRE(seams_chill_plus_graph(xyz, 1, thin, 3.5, "cutoff", &label) == 1);
  thin[4] = 10;
  thin[5] = 0;
  REQUIRE(seams_chill_plus_graph(xyz, 1, thin, 3.5, "cutoff", &label) == 1);
  REQUIRE(seams_chill_plus_rule(xyz, 1, box, 3.5, "cutoff", "nope", &label) == 1);
  REQUIRE(seams_chill_plus_rule(nullptr, 1, box, 3.5, "cutoff", "chill+", &label) == 1);
  REQUIRE(seams_chill_plus_rule(xyz, 1, box, 3.5, nullptr, "chill+", &label) == 0);
  double box9[9] = {0, 0, 0, 10, 10, 10, 0, 0, 0};
  REQUIRE(seams_chill_plus_graph_lammps(nullptr, 1, box9, 3.5, "cutoff", &label) == 1);
  REQUIRE(seams_chill_plus_graph_lammps(xyz, 1, nullptr, 3.5, "cutoff", &label) == 1);
  REQUIRE(seams_chill_plus_graph_lammps(xyz, 0, box9, 3.5, "cutoff", &label) == 1);
  REQUIRE(seams_chill_plus_graph_lammps(xyz, 1, box9, 0.0, "cutoff", &label) == 1);
  REQUIRE(seams_chill_plus_graph_lammps(xyz, 1, box9, 3.5, "nope", &label) == 1);
  REQUIRE(seams_chill_plus_graph_lammps(xyz, 1, box9, 3.5, nullptr, &label) == 0);
  double bad9[9] = {0, 0, 0, 0, 10, 10, 0, 0, 0};
  REQUIRE(seams_chill_plus_graph_lammps(xyz, 1, bad9, 3.5, "cutoff", &label) == 1);
  REQUIRE(seams_chill_plus_rule(xyz, 1, nullptr, 3.5, "cutoff", "chill+", &label) == 1);
  REQUIRE(seams_chill_plus_rule(xyz, 1, box, 3.5, "cutoff", "chill+", nullptr) == 1);
  REQUIRE(seams_chill_plus_rule(xyz, 0, box, 3.5, "cutoff", "chill+", &label) == 1);
  REQUIRE(seams_chill_plus_rule(xyz, 1, box, 0.0, "cutoff", "chill+", &label) == 1);
  REQUIRE(seams_chill_plus_rule(xyz, 1, box, 3.5, "", "chill+", &label) == 0);
  REQUIRE(seams_chill_plus_rule(xyz, 1, box, 3.5, "mutual", "chill+", &label) == 0);
  REQUIRE(seams_chill_plus_rule(xyz, 1, box, 3.5, "nope", "chill+", &label) == 1);
  thin[3] = 0;
  thin[4] = 10;
  thin[5] = 10;
  REQUIRE(seams_chill_plus_rule(xyz, 1, thin, 3.5, "cutoff", "chill+", &label) == 1);
  thin[3] = 10;
  thin[4] = 0;
  REQUIRE(seams_chill_plus_rule(xyz, 1, thin, 3.5, "cutoff", "chill+", &label) == 1);
  thin[4] = 10;
  thin[5] = 0;
  REQUIRE(seams_chill_plus_rule(xyz, 1, thin, 3.5, "cutoff", "chill+", &label) == 1);
  REQUIRE(seams_chill_plus_graph_lammps(xyz, 1, box9, 3.5, "", &label) == 0);
  REQUIRE(seams_chill_plus_graph_lammps(xyz, 1, box9, 3.5, "mutual", &label) == 0);
  REQUIRE(seams_chill_plus_graph_lammps(xyz, 1, box9, 3.5, "cutoff", nullptr) == 1);
}

TEST_CASE("CHILL and CHILL+ cover the bond, harmonic and offload arms",
          "[chill][rule]") {
  const std::array<double, 2> ang{{0.4, 1.1}};
  REQUIRE(sph::spheriHarmo(-1, ang).empty());
  REQUIRE(sph::spheriHarmo(0, ang).size() == 1);
  REQUIRE(sph::spheriHarmo(1, ang).size() == 3);
  REQUIRE(sph::spheriHarmo(2, ang).size() == 5);
  REQUIRE(sph::spheriHarmo(5, ang).size() == 11);
  REQUIRE(sph::spheriHarmo(8, ang).size() == 17);
  REQUIRE(sph::lookupTableQ3Vec(ang).size() == 7);
  REQUIRE(sph::lookupTableQ4Vec(ang).size() == 9);
  REQUIRE(sph::lookupTableQ6Vec(ang).size() == 13);
  const auto y3 = sph::lookupTableQ3Vec(ang);
  const auto y4 = sph::lookupTableQ4Vec(ang);
  const auto y6 = sph::lookupTableQ6Vec(ang);
  const auto y8 = sph::lookupTableQ8Vec(ang);
  REQUIRE(y3.size() == 7);
  REQUIRE(y4.size() == 9);
  REQUIRE(y6.size() == 13);
  REQUIRE(y8.size() == 17);
  REQUIRE(sph::lookupTableQ3(0, ang) == y3.front());
  REQUIRE(sph::lookupTableQ3(6, ang) == y3.back());
  REQUIRE(sph::lookupTableQ4(0, ang) == y4.front());
  REQUIRE(sph::lookupTableQ4(8, ang) == y4.back());
  REQUIRE(sph::lookupTableQ6(0, ang) == y6.front());
  REQUIRE(sph::lookupTableQ6(12, ang) == y6.back());
  REQUIRE(sph::lookupTableQ8(0, ang) == y8.front());
  REQUIRE(sph::lookupTableQ8(8, ang) == y8[8]);
  REQUIRE(sph::lookupTableQ8(16, ang) == y8.back());
  REQUIRE(sph::radialCoord({0.0, 0.0, 0.0})[1] == 0.0);
  const auto zenith = sph::radialCoord({0.0, 0.0, 2.0});
  REQUIRE(zenith[1] == 0.0);
  REQUIRE_THROWS_AS(chill::bondClassifier("no-such-rule"), std::out_of_range);
  std::vector<std::complex<double>> filled = sph::spheriHarmo(3, ang);
  sph::spheriHarmoInto(-1, ang, filled);
  REQUIRE(filled.empty());
  chill::registerBondClassifier("custom-chill", chill::chillPlusRule());
  REQUIRE(chill::bondClassifier("CHILL").coordinationNumber == 4);
  REQUIRE(chill::bondClassifier("custom-chill").eclipsedMax == 0.25);
  const auto axis = sph::radialCoord({1.0, 0.0, 0.0});
  REQUIRE(axis[0] != 0.0);

  Cloud hand;
  hand.nop = 6;
  hand.currentFrame = 1;
  hand.box = {30.0, 30.0, 30.0};
  hand.boxLow = {0.0, 0.0, 0.0};
  hand.pts.resize(6);
  for (int i = 0; i < hand.nop; ++i) {
    auto &p = hand.pts[static_cast<std::size_t>(i)];
    p.x = static_cast<double>(i) * 2.0;
    p.y = 0.2;
    p.z = 0.4;
    p.type = 1;
    p.atomID = i + 1;
    p.molID = p.atomID;
    p.inSlice = true;
    hand.idIndexMap[p.atomID] = i;
  }
  hand.pts[5].inSlice = false;
  auto bond = [](molSys::bond_type kind, double value) {
    molSys::Result row;
    row.classifier = kind;
    row.c_value = value;
    return row;
  };
  using Kind = molSys::bond_type;
  hand.pts[0].c_ij = {bond(Kind::staggered, -0.9), bond(Kind::staggered, -0.9),
                      bond(Kind::staggered, -0.9), bond(Kind::staggered, -0.9)};
  hand.pts[1].c_ij = {bond(Kind::staggered, -0.9), bond(Kind::staggered, -0.9),
                      bond(Kind::staggered, -0.9), bond(Kind::eclipsed, -0.1)};
  hand.pts[2].c_ij = {bond(Kind::staggered, -0.9), bond(Kind::staggered, -0.9),
                      bond(Kind::out_of_range, 0.6), bond(Kind::out_of_range, 0.7)};
  hand.pts[3].c_ij = {bond(Kind::out_of_range, 0.4), bond(Kind::out_of_range, 0.4),
                      bond(Kind::out_of_range, 0.4), bond(Kind::out_of_range, 0.4)};
  hand.pts[4].c_ij = {bond(Kind::staggered, -0.9)};
  hand.pts[5].c_ij = {bond(Kind::staggered, -0.9), bond(Kind::staggered, -0.9),
                      bond(Kind::staggered, -0.9), bond(Kind::staggered, -0.9)};
  std::vector<std::vector<int>> nList = {{1, 2, 3, 4, 5}, {2, 1, 3, 4, 5},
                                         {3, 1, 2, 4, 5}, {4, 1, 2, 3, 5},
                                         {5, 1},          {6, 1, 2, 3, 4}};
  chill::getIceTypeNoPrint(hand, nList, false);
  REQUIRE(hand.pts[0].iceType == molSys::atom_state_type::cubic);
  REQUIRE(hand.pts[1].iceType == molSys::atom_state_type::hexagonal);
  REQUIRE(hand.pts[2].iceType == molSys::atom_state_type::interfacial);
  REQUIRE(hand.pts[3].iceType == molSys::atom_state_type::water);
  REQUIRE(hand.pts[4].iceType == molSys::atom_state_type::water);
  chill::getIceTypeNoPrint(hand, nList, true);
  REQUIRE(hand.pts[5].iceType == molSys::atom_state_type::cubic);
  const auto outDir = std::filesystem::temp_directory_path() / "dseams-chill-ice";
  std::filesystem::create_directories(outDir);
  hand.pts[0].iceType = molSys::atom_state_type::cubic;
  hand.pts[1].iceType = molSys::atom_state_type::hexagonal;
  hand.pts[2].iceType = molSys::atom_state_type::interfacial;
  hand.pts[3].iceType = molSys::atom_state_type::clathrate;
  hand.pts[4].iceType = molSys::atom_state_type::interClathrate;
  hand.pts[5].iceType = molSys::atom_state_type::water;
  chill::getIceType(hand, nList, outDir.string() + "/", 1, false, "chill.txt");
  hand.pts[0].iceType = molSys::atom_state_type::cubic;
  hand.pts[1].iceType = molSys::atom_state_type::hexagonal;
  hand.pts[2].iceType = molSys::atom_state_type::interfacial;
  hand.pts[3].iceType = molSys::atom_state_type::clathrate;
  hand.pts[4].iceType = molSys::atom_state_type::interClathrate;
  hand.pts[5].iceType = molSys::atom_state_type::water;
  hand.currentFrame = 1;
  REQUIRE(chill::printIceType(hand, outDir.string() + "/", 1, false, "types.txt") >= 0);
  hand.pts[0].iceType = molSys::atom_state_type::reCubic;
  hand.pts[1].iceType = molSys::atom_state_type::reHex;
  hand.pts[2].iceType = molSys::atom_state_type::unclassified;
  hand.currentFrame = 2;
  REQUIRE(chill::printIceType(hand, outDir.string() + "/", 1, true, "types.txt") >= 0);
  std::filesystem::remove_all(outDir);

  REQUIRE(chill::isInterfacial(hand, nList, 2, 2, 0, false));
  hand.pts[2].c_ij = {bond(Kind::staggered, -0.9), bond(Kind::staggered, -0.9)};
  REQUIRE(chill::isInterfacial(hand, nList, 1, 3, 0, false));
  REQUIRE(chill::isInterfacial(hand, nList, 1, 3, 0, true));
  REQUIRE_FALSE(chill::isInterfacial(hand, nList, 0, 4, 0, true));
  REQUIRE(sph::lookupTableQ3(-1, ang) == std::complex<double>(0.0, 0.0));
  REQUIRE(sph::lookupTableQ3(99, ang) == std::complex<double>(0.0, 0.0));
  REQUIRE(sph::lookupTableQ4(99, ang) == std::complex<double>(0.0, 0.0));
  REQUIRE(sph::lookupTableQ4(-1, ang) == std::complex<double>(0.0, 0.0));
  REQUIRE(sph::lookupTableQ6(99, ang) == std::complex<double>(0.0, 0.0));
  REQUIRE(sph::lookupTableQ6(-1, ang) == std::complex<double>(0.0, 0.0));
  REQUIRE(sph::lookupTableQ8(-1, ang) == std::complex<double>(0.0, 0.0));
  REQUIRE(sph::lookupTableQ8(99, ang) == std::complex<double>(0.0, 0.0));

  std::vector<double> q6(static_cast<std::size_t>(hand.nop), 0.1);
  q6[0] = 0.8;
  q6[1] = 0.8;
  hand.pts[0].iceType = molSys::atom_state_type::water;
  hand.pts[0].c_ij = {bond(Kind::staggered, -0.9), bond(Kind::staggered, -0.9)};
  hand.pts[1].iceType = molSys::atom_state_type::water;
  hand.pts[1].c_ij = {bond(Kind::staggered, -0.8), bond(Kind::staggered, -0.8)};
  hand.pts[3].iceType = molSys::atom_state_type::cubic;
  hand.pts[4].iceType = molSys::atom_state_type::water;
  hand.pts[4].c_ij = {bond(Kind::out_of_range, -0.2), bond(Kind::out_of_range, -0.2)};
  q6[4] = 0.9;
  chill::reclassifyWater(hand, q6);
  REQUIRE(hand.pts[0].iceType == molSys::atom_state_type::reCubic);
  REQUIRE(hand.pts[1].iceType == molSys::atom_state_type::reHex);
  REQUIRE(hand.pts[3].iceType == molSys::atom_state_type::cubic);

  const auto qEmpty = chill::getq6(hand, std::vector<std::vector<int>>(6, std::vector<int>{1}), false);
  REQUIRE(qEmpty.size() == 6);
  hand.pts[0].x = hand.pts[1].x;
  hand.pts[0].y = hand.pts[1].y;
  hand.pts[0].z = hand.pts[1].z;
  const auto qSame = chill::getq6(hand, nList, true);
  REQUIRE(qSame.size() == 6);

  Cloud ice = diamond(1);
  auto rows = nneigh::neighListO(3.5, ice, 1);
  rows[0].push_back(999999);
  ice.pts[0].inSlice = false;
  struct Guard {
    seams::product::Mode mode;
    seams::product::Reduction reduction;
    Guard() : mode(seams::product::mode()), reduction(seams::product::reduction()) {}
    ~Guard() {
      seams::product::setMode(mode);
      seams::product::setReduction(reduction);
    }
  } guard;
  seams::product::setMode(seams::product::Mode::Dot2);
  seams::product::setReduction(seams::product::Reduction::Pairwise);
  chill::getCorrel(ice, rows, true, 4);
  chill::getCorrelPlus(ice, rows, true, 0, true);
  chill::getIceTypePlusNoPrint(ice, rows, true, true);
  chill::getIceType(ice, rows, outDir.string() + "/", 1, true, "slice.txt");
  seams::product::setReduction(seams::product::Reduction::Left);
  chill::getCorrel(ice, rows, false, 4);
  seams::product::setMode(seams::product::Mode::Fp64);
  seams::product::setReduction(seams::product::Reduction::Left);
  for (auto &pt : ice.pts) {
    pt.atomID += 100000;
  }
  ice.idIndexMap.clear();
  for (int i = 0; i < ice.nop; ++i) {
    ice.idIndexMap[ice.pts[static_cast<std::size_t>(i)].atomID] = i;
  }
  auto sparse = nneigh::neighListO(3.5, ice, 1);
  chill::getCorrel(ice, sparse, false, 0);
  REQUIRE(ice.pts[1].c_ij.size() > 0);
  sparse[0].push_back(1);
  chill::getCorrel(ice, sparse, false, 4);
  Cloud stuck = diamond(1);
  stuck.pts[1].x = stuck.pts[0].x;
  stuck.pts[1].y = stuck.pts[0].y;
  stuck.pts[1].z = stuck.pts[0].z;
  auto stuckRows = nneigh::neighListO(3.5, stuck, 1);
  chill::getCorrel(stuck, stuckRows, false, 4);
  Cloud dense = diamond(1);
  auto denseRows = nneigh::neighListO(3.5, dense, 1);
  dense.pts[0].atomID = -5;
  dense.idIndexMap[-5] = 0;
  denseRows[0].push_back(-1);
  denseRows[0].push_back(999999);
  chill::getCorrel(dense, denseRows, false, 4);
  chill::getCorrelPlus(dense, denseRows, false, 2);
  chill::getIceTypePlusNoPrint(dense, denseRows, false, true);

  Cloud flat;
  flat.nop = 2;
  flat.pts.resize(2);
  flat.pts[0].atomID = 1;
  flat.pts[1].atomID = 2;
  flat.pts[1].x = 1.2;
  flat.idIndexMap[1] = 0;
  flat.idIndexMap[2] = 1;
  const auto bare = chill::specializedChillPlus(flat, {});
  REQUIRE(bare.iceType.size() == 2);
  flat.box = {10.0, 10.0, 0.0};
  std::vector<std::vector<int>> pair{{1, 2}, {2, 1}};
  const auto zeroL = chill::specializedChillPlus(flat, pair);
  REQUIRE(zeroL.iceType.size() == 2);
  flat.box = {10.0};
  const auto shortBox = chill::specializedChillPlus(flat, {});
  REQUIRE(shortBox.iceType.size() == 2);
  flat.pts[1].x = 0.0;
  flat.box = {10.0, 10.0, 10.0};
  const auto samePoint = chill::specializedChillPlus(flat, pair);
  REQUIRE(samePoint.iceType.size() == 2);
  Cloud traj;
  traj = sinp::readLammpsTrj("traj/genice_sI.lammpstrj", 1, traj);
  auto trajRows = nneigh::neighListO(3.5, traj, 1);
  const auto clath = chill::specializedChillPlus(traj, trajRows);
  REQUIRE(clath.iceType.size() == static_cast<std::size_t>(traj.nop));
  Cloud mix;
  mix = sinp::readLammpsTrj("traj/exampleTraj.lammpstrj", 1, mix);
  auto mixRows = nneigh::neighListO(3.5, mix, 2);
  const auto mixed = chill::specializedChillPlus(mix, mixRows);
  REQUIRE(mixed.iceType.size() == static_cast<std::size_t>(mix.nop));
}

// A flat sheet of spacing 2.8 has opposite bonds. The l = 3 sum cancels
// and the leftover correlation is roundoff. A fractional shift must not
// move an atom across a CHILL+ class.
TEST_CASE("cancelled shells keep CHILL+ labels when the cloud is shifted", "[chill]") {
  auto sheet = [](int n) {
    Cloud cloud;
    const double spacing = 2.8;
    const double box = spacing * 8.0 + 10.0;
    cloud.box = {box, box, box};
    cloud.boxLow = {0.0, 0.0, 0.0};
    const int side = 4;
    for (int i = 0; i < n; i++) {
      molSys::Point<double> pt;
      pt.type = 1;
      pt.atomID = i + 1;
      pt.x = spacing * static_cast<double>(i % side) + 0.001 * (i + 1);
      pt.y = spacing * static_cast<double>((i / side) % side) + 0.0007 * (i + 1);
      pt.z = spacing * static_cast<double>(i / (side * side)) + 0.0005 * (i + 1);
      cloud.pts.push_back(pt);
      cloud.idIndexMap[pt.atomID] = i;
    }
    cloud.nop = n;
    return cloud;
  };
  auto labels = [](Cloud cloud) {
    const auto rows = nneigh::neighListO(3.5, cloud, 1);
    chill::applyChillPlus(cloud, rows, chill::ChillPlusRule::Published, false);
    std::vector<int> out;
    out.reserve(cloud.pts.size());
    for (const auto &pt : cloud.pts) {
      out.push_back(static_cast<int>(pt.iceType));
    }
    return out;
  };
  auto shift = [](Cloud cloud) {
    const double lx = cloud.box[0];
    const double ly = cloud.box[1];
    const double lz = cloud.box[2];
    const double add[3] = {0.02 * lx, 0.01 * ly, 0.015 * lz};
    for (auto &pt : cloud.pts) {
      double x = pt.x + add[0];
      double y = pt.y + add[1];
      double z = pt.z + add[2];
      x -= lx * std::floor(x / lx);
      y -= ly * std::floor(y / ly);
      z -= lz * std::floor(z / lz);
      pt.x = x;
      pt.y = y;
      pt.z = z;
    }
    return cloud;
  };
  for (const int n : {16, 24}) {
    const Cloud cloud = sheet(n);
    REQUIRE(labels(cloud) == labels(shift(cloud)));
  }
}

// Nguyen and Molinero, J. Phys. Chem. B 2015, 119, 9369. Equation 1
// averages q3 over the four nearest neighbours. Table 1, footnote a,
// counts the whole cutoff shell and requires four for a crystal row.
// The staggered edge is c_ij <= -0.8 and the eclipsed band is
// -0.35 <= c_ij <= 0.25. An interfacial row also looks at a neighbour:
// two staggered bonds need a neighbour with more than two, and three
// staggered bonds with no eclipsed bond need a neighbour with more
// than one. OVITO 3.17.2 averages q3 over the whole shell, uses open
// edges, and skips that neighbour test. On this oxygen the shell holds
// four neighbours, none of the four correlations sits on an edge, and
// one neighbour has a fifth oxygen inside 3.5 A. That fifth oxygen
// moves one correlation into the eclipsed band, so the published rule
// counts two eclipsed bonds (liquid) and ovito-chill+ counts three
// (interfacial clathrate).
TEST_CASE("liquid CHILL+ keeps the published four-nearest q3", "[chill][rule]") {
  auto frame = readXyz("chill-def/liquid-four-nearest.txt");
  auto rows = nneigh::neighListO(3.5, frame.cloud, 1);
  REQUIRE(chill::cutoffShellCount(frame.cloud, rows, frame.focus) == 4);
  int wider = 0;
  for (std::size_t k = 1; k < rows[static_cast<std::size_t>(frame.focus)].size(); ++k) {
    const auto found = frame.cloud.idIndexMap.find(rows[static_cast<std::size_t>(frame.focus)][k]);
    REQUIRE(found != frame.cloud.idIndexMap.end());
    if (chill::cutoffShellCount(frame.cloud, rows, found->second) != 4) {
      ++wider;
    }
  }
  REQUIRE(wider >= 1);
  auto published = frame.cloud;
  auto ovito = frame.cloud;
  chill::applyChillPlus(published, rows, chill::ChillPlusRule::Published, false);
  chill::applyChillPlus(ovito, rows, chill::ChillPlusRule::Ovito, false);
  const auto &focusP = published.pts[static_cast<std::size_t>(frame.focus)];
  const auto &focusO = ovito.pts[static_cast<std::size_t>(frame.focus)];
  REQUIRE(focusP.iceType == molSys::atom_state_type::water);
  REQUIRE(focusO.iceType == molSys::atom_state_type::interClathrate);
  auto tally = [](const std::vector<molSys::Result> &bonds) {
    int staggered = 0;
    int eclipsed = 0;
    for (const auto &bond : bonds) {
      REQUIRE(std::abs(bond.c_value + 0.8) > 1.0e-6);
      REQUIRE(std::abs(bond.c_value + 0.35) > 1.0e-6);
      REQUIRE(std::abs(bond.c_value - 0.25) > 1.0e-6);
      if (bond.classifier == molSys::bond_type::staggered) {
        ++staggered;
      } else if (bond.classifier == molSys::bond_type::eclipsed) {
        ++eclipsed;
      }
    }
    return std::pair<int, int>(staggered, eclipsed);
  };
  REQUIRE(tally(focusP.c_ij) == std::pair<int, int>(0, 2));
  REQUIRE(tally(focusO.c_ij) == std::pair<int, int>(0, 3));
  std::vector<double> left;
  std::vector<double> right;
  for (const auto &bond : focusP.c_ij) left.push_back(bond.c_value);
  for (const auto &bond : focusO.c_ij) right.push_back(bond.c_value);
  std::sort(left.begin(), left.end());
  std::sort(right.begin(), right.end());
  REQUIRE(left.size() == 4);
  REQUIRE(right.size() == 4);
  double moved = 0.0;
  for (std::size_t i = 0; i < left.size(); ++i) {
    moved = std::max(moved, std::abs(left[i] - right[i]));
  }
  REQUIRE(moved > 0.1);
}
