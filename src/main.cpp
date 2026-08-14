//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <opt_parser.h>
#include <seams_yaml.hpp>

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

#include <rang.hpp>
#include <sol/sol.hpp>

#include <iostream>
#include <string>
#include <vector>

namespace {

void registerCommon(sol::state &lua) {
  lua.set_function("readFrameOnlyOne", sinp::readLammpsTrjreduced);
  lua.set_function("readFrameOnlyOneAllAtoms", sinp::readLammpsTrj);
  lua.set_function("readFrame", sinp::readLammpsTrjO);
  lua.set_function("neighborList", nneigh::neighListO);
  lua.set_function("getHbondNetwork", bond::populateHbonds);
  lua.set_function("getHbondNetworkFromClouds",
                   bond::populateHbondsWithInputClouds);
  lua.set_function("bondNetworkByIndex", nneigh::neighbourListByIndex);
  lua.set_function("getPrimitiveRings", primitive::ringNetwork);
}

void registerDescriptors(sol::state &lua) {
  using Cloud = molSys::PointCloud<molSys::Point<double>, double>;
  lua.set_function(
      "classifyTemplates",
      [](sol::this_state ts, const Cloud &cloud,
         const std::vector<std::vector<int>> &nList, int kNeigh) {
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
  lua.set_function("soapSpectrum", chill::soapSpectrum);
  lua.set_function("soapSpectrumAll", chill::soapSpectrumAll);
  auto packQl = [](sol::this_state ts, const chill::SteinhardtQl &q) {
    sol::state_view lua(ts);
    sol::table t = lua.create_table(0, 2);
    t["ql"] = q.ql;
    t["qlBar"] = q.qlBar;
    return t;
  };
  lua.set_function("steinhardtQl",
                   [packQl](sol::this_state ts, const Cloud &cloud,
                            const std::vector<std::vector<int>> &nList,
                            int orderL) {
                     return packQl(ts, chill::steinhardtQl(cloud, nList, orderL));
                   });
  lua.set_function("steinhardtQlVoronoi",
                   [packQl](sol::this_state ts, const Cloud &cloud, double cut,
                            int orderL) {
                     return packQl(ts,
                                   chill::steinhardtQlVoronoi(cloud, cut, orderL));
                   });
  lua.set_function("voronoiFeatures", chill::voronoiFeatures);
}

} // namespace

int main(int argc, char *argv[]) {
  auto result = parse(argc, argv);
  const auto cfg = seams::loadLuaConfig(result["c"].as<std::string>());
  const std::string tFile = cfg.trajectory;

  sol::state lua;
  lua.open_libraries();

  if (cfg.topoTwoDim) {
    lua.script_file(cfg.variables);
    molSys::PointCloud<molSys::Point<double>, double> resCloud;
    std::vector<std::vector<int>> nList, hbnList;
    std::vector<std::vector<int>> rings;
    std::vector<double> rdfValues;
    auto lscript = lua.get<std::string>("functionScript");
    lua["doBOP"] = cfg.bulkUse;
    lua["topoOneDim"] = cfg.topoOneDim;
    lua["topoTwoDim"] = cfg.topoTwoDim;
    lua["topoBulk"] = cfg.bulkUse;
    lua["nList"] = &nList;
    lua["hbnList"] = &hbnList;
    lua["resCloud"] = &resCloud;
    lua["trajectory"] = tFile;
    lua["ringsAllSizes"] = &rings;
    lua["rdf"] = &rdfValues;
    registerCommon(lua);
    lua.set_function("ringAnalysis", ring::polygonRingAnalysis);
    lua.set_function("calcRDF", rdf2::rdf2Danalysis_AA);
    lua.script_file(lscript);
  }

  if (cfg.topoOneDim) {
    lua.script_file(cfg.variables);
    molSys::PointCloud<molSys::Point<double>, double> resCloud, oCloud, hCloud;
    std::vector<std::vector<int>> nList, hbnList, rings;
    int atomID = 0;
    auto lscript = lua.get<std::string>("functionScript");
    lua["doBOP"] = cfg.bulkUse;
    lua["topoOneDim"] = cfg.topoOneDim;
    lua["topoTwoDim"] = cfg.topoTwoDim;
    lua["topoBulk"] = cfg.bulkUse;
    lua["nList"] = &nList;
    lua["hbnList"] = &hbnList;
    lua["resCloud"] = &resCloud;
    lua["oCloud"] = &oCloud;
    lua["hCloud"] = &hCloud;
    lua["trajectory"] = tFile;
    lua["ringsAllSizes"] = &rings;
    lua["lowestAtomID"] = &atomID;
    registerCommon(lua);
    lua.set_function("prismAnalysis", ring::prismAnalysis);
    lua.script_file(lscript);
  }

  if (cfg.bulkUse) {
    lua.script_file(cfg.variables);
    molSys::PointCloud<molSys::Point<double>, double> resCloud, solCloud, oCloud,
        hCloud;
    std::vector<std::vector<int>> nList, hbnList, iceList, rings;
    std::vector<double> avgQ6;
    auto lscript = lua.get<std::string>("functionScript");
    lua["doBOP"] = cfg.bulkBop;
    lua["topoOneDim"] = cfg.topoOneDim;
    lua["topoTwoDim"] = cfg.topoTwoDim;
    lua["topoBulk"] = cfg.bulkTopo;
    lua["nList"] = &nList;
    lua["hbnList"] = &hbnList;
    lua["iceNeighbourList"] = &iceList;
    lua["resCloud"] = &resCloud;
    lua["oCloud"] = &oCloud;
    lua["hCloud"] = &hCloud;
    lua["clusterCloud"] = &solCloud;
    lua["avgQ6"] = &avgQ6;
    lua["trajectory"] = tFile;
    lua["ringsAllSizes"] = &rings;
    lua.set_function("writeDump", sout::writeDump);
    lua.set_function("writeHistogram", sout::writeHisto);
    registerCommon(lua);
    lua.set_function("chillPlus_cij",
                     [](molSys::PointCloud<molSys::Point<double>, double> &c,
                        const std::vector<std::vector<int>> &n, bool slice)
                         -> molSys::PointCloud<molSys::Point<double>, double> & {
                       chill::getCorrelPlus(c, n, slice);
                       return c;
                     });
    lua.set_function("chillPlus_iceType",
                     [](molSys::PointCloud<molSys::Point<double>, double> &c,
                        const std::vector<std::vector<int>> &n, std::string path,
                        int first, bool slice, std::string outName)
                         -> molSys::PointCloud<molSys::Point<double>, double> & {
                       chill::getIceTypePlus(c, n, path, first, slice, outName);
                       return c;
                     });
    lua.set_function("chill_cij",
                     [](molSys::PointCloud<molSys::Point<double>, double> &c,
                        const std::vector<std::vector<int>> &n, bool slice)
                         -> molSys::PointCloud<molSys::Point<double>, double> & {
                       chill::getCorrel(c, n, slice);
                       return c;
                     });
    lua.set_function("chill_iceType",
                     [](molSys::PointCloud<molSys::Point<double>, double> &c,
                        const std::vector<std::vector<int>> &n, std::string path,
                        int first, bool slice, std::string outName)
                         -> molSys::PointCloud<molSys::Point<double>, double> & {
                       chill::getIceType(c, n, path, first, slice, outName);
                       return c;
                     });
    lua.set_function("averageQ6", chill::getq6);
    lua.set_function("modifyChill",
                     [](molSys::PointCloud<molSys::Point<double>, double> &c,
                        std::vector<double> &q6)
                         -> molSys::PointCloud<molSys::Point<double>, double> & {
                       chill::reclassifyWater(c, q6);
                       return c;
                     });
    lua.set_function("percentage_Ice", chill::printIceType);
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
    lua.set_function("bulkTopologicalNetworkCriterion",
                     ring::topoBulkAnalysis);
    lua.set_function("bulkTopoUnitMatching", tum3::topoUnitMatchingBulk);
    lua.script_file(lscript);
  }

  if (cfg.structureDesc) {
    lua.script_file(cfg.variables);
    molSys::PointCloud<molSys::Point<double>, double> resCloud;
    std::vector<std::vector<int>> nList;
    auto lscript = lua.get<std::string>("functionScript");
    lua["resCloud"] = &resCloud;
    lua["nList"] = &nList;
    lua["trajectory"] = tFile;
    registerCommon(lua);
    registerDescriptors(lua);
    lua.script_file(lscript);
  }

  std::cout << rang::style::bold << "Welcome to the Black Parade.\nYou ran:-\n"
            << rang::style::reset
            << "\nBulk Ice Analysis: " << (cfg.bulkUse ? "true" : "false")
            << "\nQuasi-one-dimensional Ice Analysis: "
            << (cfg.topoOneDim ? "true" : "false")
            << "\nQuasi-two-dimensional Ice Analysis: "
            << (cfg.topoTwoDim ? "true" : "false")
            << "\nStructure descriptors: "
            << (cfg.structureDesc ? "true" : "false") << "\n";
  return 0;
}
