//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <lua_api.hpp>

#include <bond.hpp>
#include <bop.hpp>
#include <bulkTUM.hpp>
#include <cluster.hpp>
#include <franzblau.hpp>
#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <rdf2d.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>
#include <selection.hpp>
#include <structure_desc.hpp>
#include <topo_bulk.hpp>
#include <topo_one_dim.hpp>
#include <topo_two_dim.hpp>
#include <voronoi_qlm.hpp>

#include <array>
#include <string>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

constexpr std::array<double, 3> zeroBounds{0.0, 0.0, 0.0};

//! Human-readable name of a per-particle CHILL/CHILL+ state
const char *iceStateName(molSys::atom_state_type state) {
  switch (state) {
  case molSys::atom_state_type::cubic:
    return "cubic";
  case molSys::atom_state_type::hexagonal:
    return "hexagonal";
  case molSys::atom_state_type::water:
    return "water";
  case molSys::atom_state_type::interfacial:
    return "interfacial";
  case molSys::atom_state_type::clathrate:
    return "clathrate";
  case molSys::atom_state_type::interClathrate:
    return "interClathrate";
  case molSys::atom_state_type::reCubic:
    return "reCubic";
  case molSys::atom_state_type::reHex:
    return "reHex";
  case molSys::atom_state_type::unclassified:
    break;
  }
  return "unclassified";
}

//! Per-particle ice states of a cloud, as strings ready for a Lua table
std::vector<std::string> iceStateNames(const Cloud &yCloud) {
  std::vector<std::string> names;
  names.reserve(yCloud.pts.size());
  for (const auto &pt : yCloud.pts) {
    names.emplace_back(iceStateName(pt.iceType));
  }
  return names;
}

//! Packs a SteinhardtQl result into a {ql = ..., qlBar = ...} Lua table
sol::table packSteinhardt(sol::state_view lua, const chill::SteinhardtQl &q) {
  sol::table t = lua.create_table(0, 2);
  t["ql"] = sol::as_table(q.ql);
  t["qlBar"] = sol::as_table(q.qlBar);
  return t;
}

} // namespace

namespace luaApi {

void registerIO(sol::state &lua) {
  lua.new_usertype<Cloud>(
      "PointCloud", sol::constructors<Cloud()>(), "nop",
      sol::readonly(&Cloud::nop), "currentFrame",
      sol::readonly(&Cloud::currentFrame), "box",
      [](const Cloud &c) { return sol::as_table(c.box); }, "boxLow",
      [](const Cloud &c) { return sol::as_table(c.boxLow); }, "iceTypes",
      [](const Cloud &c) { return sol::as_table(iceStateNames(c)); });
  lua.set_function(
      "readLammpsTrjO",
      [](std::string filename, int targetFrame, int typeO,
         sol::optional<bool> isSlice, sol::optional<std::array<double, 3>> low,
         sol::optional<std::array<double, 3>> high) {
        Cloud scratch;
        return sinp::readLammpsTrjO(filename, targetFrame, scratch, typeO,
                                    isSlice.value_or(false),
                                    low.value_or(zeroBounds),
                                    high.value_or(zeroBounds));
      });
  lua.set_function(
      "readLammpsTrjreduced",
      [](std::string filename, int targetFrame, int typeI,
         sol::optional<bool> isSlice, sol::optional<std::array<double, 3>> low,
         sol::optional<std::array<double, 3>> high) {
        Cloud scratch;
        return sinp::readLammpsTrjreduced(filename, targetFrame, scratch, typeI,
                                          isSlice.value_or(false),
                                          low.value_or(zeroBounds),
                                          high.value_or(zeroBounds));
      });
  lua.set_function("readXYZ", sinp::readXYZ);
#ifdef SEAMS_HAS_CHEMFILES
  lua.set_function("readChemfiles",
                   [](std::string filename, int targetFrame,
                      sol::optional<int> typeFilter) {
                     Cloud scratch;
                     return sinp::readChemfiles(filename, targetFrame, scratch,
                                                typeFilter.value_or(-1));
                   });
#endif
#ifdef SEAMS_HAS_READCON
  lua.set_function("readCon", [](std::string filename, int targetFrame) {
    Cloud scratch;
    return sinp::readCon(filename, targetFrame, scratch);
  });
#endif
  // Legacy spellings, container-userdata semantics
  lua.set_function("readFrameOnlyOne", sinp::readLammpsTrjreduced);
  lua.set_function("readFrameOnlyOneAllAtoms", sinp::readLammpsTrj);
  lua.set_function("readFrame", sinp::readLammpsTrjO);
  lua.set_function("writeDump", sout::writeDump);
  lua.set_function("writeHistogram", sout::writeHisto);
}

void registerNeighbours(sol::state &lua) {
  lua.set_function("neighListO",
                   [](double rcutoff, const Cloud &yCloud, int typeI) {
                     return sol::as_nested(
                         nneigh::neighListO(rcutoff, yCloud, typeI));
                   });
  lua.set_function("getNewNeighbourListByIndex",
                   [](const Cloud &yCloud, double cutoff) {
                     return sol::as_nested(
                         nneigh::getNewNeighbourListByIndex(yCloud, cutoff));
                   });
  // New-style bindings take the neighbour list by value: sol2 can build a
  // container from a plain Lua table only when the parameter owns it, while
  // a reference parameter binds container userdata alone.
  lua.set_function(
      "neighbourListByIndex",
      [](const Cloud &yCloud, std::vector<std::vector<int>> nList) {
        return sol::as_nested(nneigh::neighbourListByIndex(yCloud, nList));
      });
  // Legacy spellings, container-userdata semantics
  lua.set_function("neighborList", nneigh::neighListO);
  lua.set_function("bondNetworkByIndex", nneigh::neighbourListByIndex);
  // distCutoff/angleCutoff default to the water criterion (2.42 A, 30 deg)
  // Neighbour lists arrive as plain Lua tables; container parameters must be
  // taken by value or sol2 dereferences a null userdata
  lua.set_function("getHbondNetwork",
                   [](std::string filename, Cloud &yCloud,
                      std::vector<std::vector<int>> nList,
                      int targetFrame, int Htype, sol::optional<double> dist,
                      sol::optional<double> angle) {
                     return bond::populateHbonds(filename, yCloud, nList,
                                                 targetFrame, Htype,
                                                 dist.value_or(2.42),
                                                 angle.value_or(30.0));
                   });
  lua.set_function("getHbondNetworkFromClouds",
                   [](Cloud &yCloud, Cloud &hCloud,
                      std::vector<std::vector<int>> nList,
                      sol::optional<double> dist, sol::optional<double> angle) {
                     return bond::populateHbondsWithInputClouds(
                         yCloud, hCloud, nList, dist.value_or(2.42),
                         angle.value_or(30.0));
                   });
}

void registerRings(sol::state &lua) {
  lua.set_function("ringNetwork",
                   [](std::vector<std::vector<int>> nList, int maxDepth) {
                     return sol::as_nested(primitive::ringNetwork(nList, maxDepth));
                   });
  lua.new_usertype<primitive::RingUpdater>(
      "RingUpdater", sol::constructors<primitive::RingUpdater(int)>(), "update",
      [](primitive::RingUpdater &self, std::vector<std::vector<int>> nList) {
        return sol::as_nested(self.update(nList));
      },
      "lastRecomputedSources", &primitive::RingUpdater::lastRecomputedSources,
      "lastBallsRefreshed", &primitive::RingUpdater::lastBallsRefreshed);
  // Legacy spelling, container-userdata semantics
  lua.set_function("getPrimitiveRings", primitive::ringNetwork);
}

void registerOrder(sol::state &lua) {
  // Bond-classification rule sets: CHILL and CHILL+ are registered
  // instances; scripts register their own materials by name
  const auto ruleFromTable = [](const sol::table &t) {
    chill::BondClassifier rule = chill::chillRule();
    rule.staggeredMax = t.get_or("staggeredMax", rule.staggeredMax);
    rule.eclipsedMin = t.get_or("eclipsedMin", rule.eclipsedMin);
    rule.eclipsedMax = t.get_or("eclipsedMax", rule.eclipsedMax);
    rule.coordinationNumber =
        t.get_or("coordinationNumber", rule.coordinationNumber);
    return rule;
  };
  lua.set_function("classifyBonds",
                   [ruleFromTable](Cloud &yCloud,
                                   std::vector<std::vector<int>> nList,
                                   const sol::object &ruleSpec,
                                   sol::optional<bool> isSlice) {
                     const chill::BondClassifier rule =
                         ruleSpec.is<std::string>()
                             ? chill::bondClassifier(
                                   ruleSpec.as<std::string>())
                             : ruleFromTable(ruleSpec.as<sol::table>());
                     chill::classifyBonds(yCloud, nList, rule,
                                          isSlice.value_or(false));
                   });
  lua.set_function("registerBondClassifier",
                   [ruleFromTable](std::string name, const sol::table &t) {
                     chill::registerBondClassifier(name, ruleFromTable(t));
                   });
  lua.set_function("bondClassifierNames", []() {
    return sol::as_table(chill::bondClassifierNames());
  });
  lua.set_function("getCorrelPlus",
                   [](Cloud &yCloud, std::vector<std::vector<int>> nList,
                      sol::optional<bool> isSlice,
                      sol::optional<int> coordinationNumber) {
                     chill::getCorrelPlus(yCloud, nList,
                                          isSlice.value_or(false),
                                          coordinationNumber.value_or(4));
                   });
  lua.set_function(
      "getIceTypePlus",
      [](Cloud &yCloud, std::vector<std::vector<int>> nList, std::string path,
         int firstFrame, sol::optional<bool> isSlice,
         sol::optional<std::string> outputFileName) {
        chill::getIceTypePlus(yCloud, nList, path, firstFrame,
                              isSlice.value_or(false),
                              outputFileName.value_or("chillPlus.txt"));
        return sol::as_table(iceStateNames(yCloud));
      });
  lua.set_function("getCorrel",
                   [](Cloud &yCloud, std::vector<std::vector<int>> nList,
                      sol::optional<bool> isSlice,
                      sol::optional<int> coordinationNumber) {
                     chill::getCorrel(yCloud, nList, isSlice.value_or(false),
                                      coordinationNumber.value_or(4));
                   });
  lua.set_function(
      "getIceType",
      [](Cloud &yCloud, std::vector<std::vector<int>> nList, std::string path,
         int firstFrame, sol::optional<bool> isSlice,
         sol::optional<std::string> outputFileName) {
        chill::getIceType(yCloud, nList, path, firstFrame,
                          isSlice.value_or(false),
                          outputFileName.value_or("chill.txt"));
        return sol::as_table(iceStateNames(yCloud));
      });
  lua.set_function("steinhardtQl",
                   [](sol::this_state ts, const Cloud &yCloud,
                      std::vector<std::vector<int>> nList, int orderL) {
                     return packSteinhardt(
                         sol::state_view(ts),
                         chill::steinhardtQl(yCloud, nList, orderL));
                   });
  lua.set_function("steinhardtQlVoronoi",
                   [](sol::this_state ts, const Cloud &yCloud,
                      double candidateCutoff, int orderL) {
                     return packSteinhardt(sol::state_view(ts),
                                           chill::steinhardtQlVoronoi(
                                               yCloud, candidateCutoff, orderL));
                   });
  lua.set_function(
      "voronoiFacetWeights",
      [](sol::this_state ts, const Cloud &yCloud, double candidateCutoff) {
        sol::state_view lua(ts);
        const auto cells = chill::voronoiFacetWeights(yCloud, candidateCutoff);
        sol::table out = lua.create_table(static_cast<int>(cells.size()), 0);
        for (size_t i = 0; i < cells.size(); i++) {
          sol::table cell = lua.create_table(0, 2);
          cell["neighbours"] = sol::as_table(cells[i].neighbours);
          cell["weights"] = sol::as_table(cells[i].weights);
          out[i + 1] = cell;
        }
        return out;
      });
  // Legacy spellings kept for the bulk example scripts
  lua.set_function("chillPlus_cij",
                   [](Cloud &c, const std::vector<std::vector<int>> &n,
                      bool slice) -> Cloud & {
                     chill::getCorrelPlus(c, n, slice);
                     return c;
                   });
  lua.set_function("chillPlus_iceType",
                   [](Cloud &c, const std::vector<std::vector<int>> &n,
                      std::string path, int first, bool slice,
                      std::string outName) -> Cloud & {
                     chill::getIceTypePlus(c, n, path, first, slice, outName);
                     return c;
                   });
  lua.set_function("chill_cij",
                   [](Cloud &c, const std::vector<std::vector<int>> &n,
                      bool slice) -> Cloud & {
                     chill::getCorrel(c, n, slice);
                     return c;
                   });
  lua.set_function("chill_iceType",
                   [](Cloud &c, const std::vector<std::vector<int>> &n,
                      std::string path, int first, bool slice,
                      std::string outName) -> Cloud & {
                     chill::getIceType(c, n, path, first, slice, outName);
                     return c;
                   });
  lua.set_function("averageQ6", chill::getq6);
  lua.set_function("modifyChill",
                   [](Cloud &c, std::vector<double> &q6) -> Cloud & {
                     chill::reclassifyWater(c, q6);
                     return c;
                   });
  lua.set_function("percentage_Ice", chill::printIceType);
}

void registerDescriptors(sol::state &lua) {
  lua.set_function(
      "classifyTemplates",
      [](sol::this_state ts, const Cloud &cloud,
         std::vector<std::vector<int>> nList, int kNeigh) {
        sol::state_view lua(ts);
        const auto hits = chill::classifyTemplates(cloud, nList, kNeigh);
        sol::table out = lua.create_table(static_cast<int>(hits.size()), 0);
        for (size_t i = 0; i < hits.size(); i++) {
          sol::table row = lua.create_table(0, 2);
          row["name"] = hits[i].name;
          row["rmsd"] = hits[i].rmsd;
          out[i + 1] = row;
        }
        return out;
      });
  lua.set_function("soapSpectrum",
                   [](const Cloud &yCloud, int iatom,
                      std::vector<std::vector<int>> nList, int nMax, int lMax,
                      double rcut) {
                     return sol::as_table(chill::soapSpectrum(
                         yCloud, iatom, nList, nMax, lMax, rcut));
                   });
  lua.set_function("soapSpectrumAll",
                   [](const Cloud &yCloud,
                      std::vector<std::vector<int>> nList, int nMax,
                      int lMax, double rcut) {
                     return sol::as_nested(chill::soapSpectrumAll(yCloud, nList,
                                                               nMax, lMax,
                                                               rcut));
                   });
  lua.set_function("voronoiFeatures",
                   [](const Cloud &yCloud, double candidateCutoff) {
                     return sol::as_nested(
                         chill::voronoiFeatures(yCloud, candidateCutoff));
                   });
}

void registerTopology(sol::state &lua) {
  lua.set_function("ringAnalysis", ring::polygonRingAnalysis);
  lua.set_function("calcRDF", rdf2::rdf2Danalysis_AA);
  lua.set_function(
      "prismAnalysis",
      [](std::string path, const std::vector<std::vector<int>> &rings,
         const std::vector<std::vector<int>> &nList, Cloud &cloud, int maxDepth,
         int atomID, int firstFrame, int currentFrame, bool doShapeMatching) {
        return ring::prismAnalysis(path, rings, nList, cloud, maxDepth, atomID,
                                   firstFrame, currentFrame, doShapeMatching);
      });
  lua.set_function("clusterAnalysis", clump::clusterAnalysis);
  lua.set_function("recenterCluster", clump::recenterClusterCloud);
  lua.set_function("getPointCloudAtomsOfOneAtomType",
                   gen::getPointCloudOneAtomType);
  lua.set_function("selectInSingleSlice", gen::moleculesInSingleSlice);
  lua.set_function("selectEdgeAtomsInRingsWithinSlice",
                   ring::getEdgeMoleculesInRings);
  lua.set_function("selectAtomsInSliceWithRingEdgeAtoms",
                   ring::printSliceGetEdgeMoleculesInRings);
  lua.set_function("bulkRingNumberAnalysis", ring::bulkPolygonRingAnalysis);
  lua.set_function("bulkTopologicalNetworkCriterion", ring::topoBulkAnalysis);
  lua.set_function("bulkTopoUnitMatching", tum3::topoUnitMatchingBulk);
}

void registerAll(sol::state &lua) {
  registerIO(lua);
  registerNeighbours(lua);
  registerRings(lua);
  registerOrder(lua);
  registerDescriptors(lua);
  registerTopology(lua);
}

} // namespace luaApi
