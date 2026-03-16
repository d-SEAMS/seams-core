/*
** This file is part of d-SEAMS (PydSEAMSlib).
**
** SPDX-License-Identifier: MIT
**
** Copyright (c) 2023--present, d-SEAMS core team
** All rights reserved.
*/
#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/unordered_map.h>

#include <bond.hpp>
#include <bop.hpp>
#include <bulkTUM.hpp>
#include <cluster.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <rdf2d.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>
#include <selection.hpp>
#include <topo_one_dim.hpp>
#include <topo_two_dim.hpp>

#include <cstdint>
#include <format>

namespace nb = nanobind;

NB_MODULE(_core, m) {
    m.doc() = "PydSEAMSlib nanobind bindings";

    nb::class_<molSys::Point<double>>(m, "PointDouble")
        .def(nb::init<>())
        .def_rw("c_type", &molSys::Point<double>::type)
        .def_rw("molID", &molSys::Point<double>::molID)
        .def_rw("atomID", &molSys::Point<double>::atomID)
        .def_rw("iceType", &molSys::Point<double>::iceType)
        .def_rw("x", &molSys::Point<double>::x)
        .def_rw("y", &molSys::Point<double>::y)
        .def_rw("z", &molSys::Point<double>::z)
        .def_rw("c_ij", &molSys::Point<double>::c_ij)
        .def_rw("inSlice", &molSys::Point<double>::inSlice)
        .def("__repr__",
             [](const molSys::Point<double> &self_C) {
                 std::uintptr_t ptr_val = std::uintptr_t(&self_C);
                 return std::format("<PointDouble mem_loc:{:x}>", static_cast<uint64_t>(ptr_val));
             })
        .def("__str__", [](const molSys::Point<double> &self_C) {
            return std::format("x: {} y: {} z: {} type: {} molID: {} atomID: {} inSlice: {}",
                               self_C.x, self_C.y, self_C.z,
                               self_C.type, self_C.molID, self_C.atomID, self_C.inSlice);
        });

    nb::enum_<molSys::bond_type>(m, "BondType")
        .value("staggered", molSys::bond_type::staggered)
        .value("eclipsed", molSys::bond_type::eclipsed)
        .value("out_of_range", molSys::bond_type::out_of_range);

    nb::enum_<molSys::atom_state_type>(m, "AtomStateType")
        .value("cubic", molSys::atom_state_type::cubic)
        .value("hexagonal", molSys::atom_state_type::hexagonal)
        .value("water", molSys::atom_state_type::water)
        .value("interfacial", molSys::atom_state_type::interfacial)
        .value("clathrate", molSys::atom_state_type::clathrate)
        .value("interClathrate", molSys::atom_state_type::interClathrate)
        .value("unclassified", molSys::atom_state_type::unclassified)
        .value("reCubic", molSys::atom_state_type::reCubic)
        .value("reHex", molSys::atom_state_type::reHex);

    nb::class_<molSys::Result>(m, "Result")
        .def(nb::init<>())
        .def_rw("classifier", &molSys::Result::classifier)
        .def_rw("c_value", &molSys::Result::c_value)
        .def("__repr__", [](const molSys::Result &self_C) {
            std::uintptr_t ptr_val = std::uintptr_t(&self_C);
            return std::format("<Result mem_loc:{:x}>", static_cast<uint64_t>(ptr_val));
        });

    nb::class_<molSys::PointCloud<molSys::Point<double>, double>>(m, "PointCloudDouble")
        .def(nb::init<>())
        .def_rw("pts", &molSys::PointCloud<molSys::Point<double>, double>::pts)
        .def_rw("currentFrame", &molSys::PointCloud<molSys::Point<double>, double>::currentFrame)
        .def_rw("nop", &molSys::PointCloud<molSys::Point<double>, double>::nop)
        .def_rw("box", &molSys::PointCloud<molSys::Point<double>, double>::box)
        .def_rw("boxLow", &molSys::PointCloud<molSys::Point<double>, double>::boxLow)
        .def_rw("idIndexMap", &molSys::PointCloud<molSys::Point<double>, double>::idIndexMap);

    // I/O (lambdas hide the yCloud* in/out parameter, creating it internally)
    m.def("readXYZ", &sinp::readXYZ, nb::arg("filename"));
    m.def("readLammpsTrjreduced",
          [](std::string filename, int targetFrame, int typeI, bool isSlice,
             std::array<double, 3> coordLow, std::array<double, 3> coordHigh) {
            molSys::PointCloud<molSys::Point<double>, double> yCloud;
            return sinp::readLammpsTrjreduced(filename, targetFrame, &yCloud,
                                               typeI, isSlice, coordLow, coordHigh);
          },
          nb::arg("filename"), nb::arg("targetFrame"), nb::arg("typeI"),
          nb::arg("isSlice"), nb::arg("coordLow"), nb::arg("coordHigh"));
    m.def("readLammpsTrjO",
          [](std::string filename, int targetFrame, int typeO, bool isSlice,
             std::array<double, 3> coordLow, std::array<double, 3> coordHigh) {
            molSys::PointCloud<molSys::Point<double>, double> yCloud;
            return sinp::readLammpsTrjO(filename, targetFrame, &yCloud,
                                         typeO, isSlice, coordLow, coordHigh);
          },
          nb::arg("filename"), nb::arg("targetFrame"), nb::arg("typeO"),
          nb::arg("isSlice"), nb::arg("coordLow"), nb::arg("coordHigh"));
    m.def("readLammpsTrj",
          [](std::string filename, int targetFrame, bool isSlice,
             std::array<double, 3> coordLow, std::array<double, 3> coordHigh) {
            molSys::PointCloud<molSys::Point<double>, double> yCloud;
            return sinp::readLammpsTrj(filename, targetFrame, &yCloud,
                                        isSlice, coordLow, coordHigh);
          },
          nb::arg("filename"), nb::arg("targetFrame"),
          nb::arg("isSlice"), nb::arg("coordLow"), nb::arg("coordHigh"));
    m.def("readBonds", &sinp::readBonds, nb::arg("filename"));
    m.def("atomInSlice", &sinp::atomInSlice,
          nb::arg("x"), nb::arg("y"), nb::arg("z"),
          nb::arg("coordLow"), nb::arg("coordHigh"));

    // Neighbours
    m.def("clearNeighbourList", &nneigh::clearNeighbourList, nb::arg("nList"));
    m.def("getNewNeighbourListByIndex", &nneigh::getNewNeighbourListByIndex,
          nb::arg("yCloud"), nb::arg("cutoff"));
    m.def("halfNeighList", &nneigh::halfNeighList,
          nb::arg("yCloud"), nb::arg("rcutoff"), nb::arg("typeI"));
    m.def("neighbourListByIndex", &nneigh::neighbourListByIndex,
          nb::arg("yCloud"), nb::arg("nList"));
    m.def("neighList", &nneigh::neighList,
          nb::arg("yCloud"), nb::arg("rcutoff"), nb::arg("typeI"), nb::arg("typeJ"));
    m.def("neighListO", &nneigh::neighListO,
          nb::arg("rcutoff"), nb::arg("yCloud"), nb::arg("typeI"));

    // Bonds
    m.def("createBondsFromCages", &bond::createBondsFromCages,
          nb::arg("rings"), nb::arg("cageList"), nb::arg("type"), nb::arg("nRings"));
    m.def("getHbondDistanceOH", &bond::getHbondDistanceOH,
          nb::arg("oCloud"), nb::arg("hCloud"), nb::arg("oAtomIndex"), nb::arg("hAtomIndex"));
    m.def("populateHbonds", &bond::populateHbonds,
          nb::arg("filename"), nb::arg("yCloud"), nb::arg("nList"),
          nb::arg("targetFrame"), nb::arg("Htype"));
    m.def("populateHbondsWithInputClouds", &bond::populateHbondsWithInputClouds,
          nb::arg("yCloud"), nb::arg("hCloud"), nb::arg("nList"));
    m.def("trimBonds", &bond::trimBonds, nb::arg("bonds"));

    // Ring analysis (Franzblau)
    m.def("clearGraph", &primitive::clearGraph, nb::arg("currentGraph"));
    m.def("countAllRingsFromIndex", &primitive::countAllRingsFromIndex,
          nb::arg("neighHbondList"), nb::arg("maxDepth"));
    m.def("ringNetwork", &primitive::ringNetwork,
          nb::arg("nList"), nb::arg("maxDepth"));
    m.def("populateGraphFromIndices", &primitive::populateGraphFromIndices, nb::arg("nList"));
    m.def("populateGraphFromNListID", &primitive::populateGraphFromNListID,
          nb::arg("yCloud"), nb::arg("neighHbondList"));
    m.def("removeNonSPrings", &primitive::removeNonSPrings, nb::arg("fullGraph"));
    m.def("restoreEdgesFromIndices", &primitive::restoreEdgesFromIndices,
          nb::arg("fullGraph"), nb::arg("nList"));

    // Ring classification
    m.def("assignPolygonType", &ring::assignPolygonType,
          nb::arg("rings"), nb::arg("atomTypes"), nb::arg("nRings"));
    m.def("assignPrismType", &ring::assignPrismType,
          nb::arg("rings"), nb::arg("listPrism"), nb::arg("ringSize"),
          nb::arg("ringType"), nb::arg("atomTypes"), nb::arg("atomState"));
    m.def("clearRingList", &ring::clearRingList, nb::arg("rings"));
    m.def("compareRings", &ring::compareRings, nb::arg("ring1"), nb::arg("ring2"));
    m.def("commonElementsInThreeRings", &ring::commonElementsInThreeRings,
          nb::arg("ring1"), nb::arg("ring2"), nb::arg("ring3"));
    m.def("deformedPrismTypes", &ring::deformedPrismTypes,
          nb::arg("atomState"), nb::arg("atomTypes"), nb::arg("maxDepth"));
    m.def("discardExtraTetragonBlocks", &ring::discardExtraTetragonBlocks,
          nb::arg("basal1"), nb::arg("basal2"), nb::arg("yCloud"));
    m.def("findPrisms", &ring::findPrisms,
          nb::arg("rings"), nb::arg("ringType"), nb::arg("nPerfectPrisms"),
          nb::arg("nImperfectPrisms"), nb::arg("nList"), nb::arg("rmsdPerAtom"),
          nb::arg("doShapeMatching"), nb::arg("yCloud"));
    m.def("findsCommonElements", &ring::findsCommonElements,
          nb::arg("ring1"), nb::arg("ring2"));
    m.def("findTripletInRing", &ring::findTripletInRing,
          nb::arg("ring"), nb::arg("triplet"));
    m.def("getSingleRingSize", &ring::getSingleRingSize,
          nb::arg("rings"), nb::arg("ringSize"));
    m.def("hasCommonElements", &ring::hasCommonElements,
          nb::arg("ring1"), nb::arg("ring2"));
    m.def("basalPrismConditions", &ring::basalPrismConditions,
          nb::arg("nList"), nb::arg("basal1"), nb::arg("basal2"));
    m.def("relaxedPrismConditions", &ring::relaxedPrismConditions,
          nb::arg("nList"), nb::arg("basal1"), nb::arg("basal2"));

    // Topology analysis
    m.def("polygonRingAnalysis", &ring::polygonRingAnalysis,
          nb::arg("path"), nb::arg("rings"), nb::arg("nList"), nb::arg("yCloud"),
          nb::arg("maxDepth"), nb::arg("sheetArea"), nb::arg("firstFrame"));
    m.def("bulkPolygonRingAnalysis", &ring::bulkPolygonRingAnalysis,
          nb::arg("path"), nb::arg("rings"), nb::arg("nList"), nb::arg("yCloud"),
          nb::arg("maxDepth"), nb::arg("firstFrame"));
    m.def("prismAnalysis", &ring::prismAnalysis,
          nb::arg("path"), nb::arg("rings"), nb::arg("nList"), nb::arg("yCloud"),
          nb::arg("maxDepth"), nb::arg("atomID"), nb::arg("firstFrame"),
          nb::arg("currentFrame"), nb::arg("doShapeMatching"));
    m.def("rmAxialTranslations", &ring::rmAxialTranslations,
          nb::arg("yCloud"), nb::arg("atomID"), nb::arg("firstFrame"), nb::arg("currentFrame"));

    // TUM (Topological Unit Matching)
    m.def("atomsFromCages", &tum3::atomsFromCages,
          nb::arg("rings"), nb::arg("cageList"), nb::arg("clusterCages"));
    m.def("averageRMSDatom", &tum3::averageRMSDatom,
          nb::arg("rmsdPerAtom"), nb::arg("noOfCommonAtoms"));
    m.def("buildRefDDC", &tum3::buildRefDDC, nb::arg("fileName"));
    m.def("buildRefHC", &tum3::buildRefHC, nb::arg("fileName"));
    m.def("clusterCages", &tum3::clusterCages,
          nb::arg("yCloud"), nb::arg("path"), nb::arg("rings"), nb::arg("cageList"),
          nb::arg("numHC"), nb::arg("numDDC"));
    m.def("shapeMatchDDC", &tum3::shapeMatchDDC,
          nb::arg("yCloud"), nb::arg("refPoints"), nb::arg("cageList"),
          nb::arg("cageIndex"), nb::arg("rings"), nb::arg("quat"), nb::arg("rmsd"));
    m.def("shapeMatchHC", &tum3::shapeMatchHC,
          nb::arg("yCloud"), nb::arg("refPoints"), nb::arg("cageUnit"),
          nb::arg("rings"), nb::arg("nList"), nb::arg("quat"), nb::arg("rmsd"));
    m.def("topoBulkCriteria", &tum3::topoBulkCriteria,
          nb::arg("path"), nb::arg("rings"), nb::arg("nList"), nb::arg("yCloud"),
          nb::arg("firstFrame"), nb::arg("numHC"), nb::arg("numDDC"), nb::arg("ringType"));
    m.def("topoUnitMatchingBulk", &tum3::topoUnitMatchingBulk,
          nb::arg("path"), nb::arg("rings"), nb::arg("nList"), nb::arg("yCloud"),
          nb::arg("firstFrame"), nb::arg("printClusters"), nb::arg("onlyTetrahedral"));
    m.def("updateRMSDatom", &tum3::updateRMSDatom,
          nb::arg("rings"), nb::arg("cageUnit"), nb::arg("rmsd"),
          nb::arg("rmsdPerAtom"), nb::arg("noOfCommonAtoms"), nb::arg("atomTypes"));

    // Selection
    m.def("getPointCloudOneAtomType", &gen::getPointCloudOneAtomType,
          nb::arg("yCloud"), nb::arg("outCloud"), nb::arg("atomTypeI"),
          nb::arg("isSlice"), nb::arg("coordLow"), nb::arg("coordHigh"));
    m.def("atomsInSingleSlice", &gen::atomsInSingleSlice,
          nb::arg("yCloud"), nb::arg("clearPreviousSliceSelection"),
          nb::arg("coordLow"), nb::arg("coordHigh"));
    m.def("moleculesInSingleSlice", &gen::moleculesInSingleSlice,
          nb::arg("yCloud"), nb::arg("clearPreviousSliceSelection"),
          nb::arg("coordLow"), nb::arg("coordHigh"));
    m.def("setAtomsWithSameMolID", &gen::setAtomsWithSameMolID,
          nb::arg("yCloud"), nb::arg("molIDAtomIDmap"), nb::arg("molID"), nb::arg("inSliceValue"));
    m.def("getEdgeMoleculesInRings", &ring::getEdgeMoleculesInRings,
          nb::arg("rings"), nb::arg("oCloud"), nb::arg("yCloud"),
          nb::arg("identicalCloud"), nb::arg("coordLow"), nb::arg("coordHigh"));
    m.def("printSliceGetEdgeMoleculesInRings", &ring::printSliceGetEdgeMoleculesInRings,
          nb::arg("path"), nb::arg("rings"), nb::arg("oCloud"), nb::arg("yCloud"),
          nb::arg("coordLow"), nb::arg("coordHigh"), nb::arg("identicalCloud"));

    // CHILL/CHILL+ classification
    m.def("getCorrelPlus", &chill::getCorrelPlus,
          nb::arg("yCloud"), nb::arg("nList"), nb::arg("isSlice"));
    m.def("getIceTypePlus", &chill::getIceTypePlus,
          nb::arg("yCloud"), nb::arg("nList"), nb::arg("path"),
          nb::arg("firstFrame"), nb::arg("isSlice"), nb::arg("outputFileName"));
    m.def("getq6", &chill::getq6,
          nb::arg("yCloud"), nb::arg("nList"), nb::arg("isSlice"));
    m.def("reclassifyWater", &chill::reclassifyWater,
          nb::arg("yCloud"), nb::arg("q6"));
    m.def("printIceType", &chill::printIceType,
          nb::arg("yCloud"), nb::arg("path"), nb::arg("firstFrame"),
          nb::arg("isSlice"), nb::arg("outputFileName"));

    // Output
    m.def("writeDump", &sout::writeDump,
          nb::arg("yCloud"), nb::arg("path"), nb::arg("outFile"));

    // Clustering
    m.def("clusterAnalysis", &clump::clusterAnalysis,
          nb::arg("path"), nb::arg("iceCloud"), nb::arg("yCloud"), nb::arg("nList"),
          nb::arg("iceNeighbourList"), nb::arg("cutoff"), nb::arg("firstFrame"),
          nb::arg("bopAnalysis"));
    m.def("recenterClusterCloud", &clump::recenterClusterCloud,
          nb::arg("iceCloud"), nb::arg("nList"));

    // RDF
    m.def("rdf2Danalysis_AA", &rdf2::rdf2Danalysis_AA,
          nb::arg("path"), nb::arg("rdfValues"), nb::arg("yCloud"),
          nb::arg("cutoff"), nb::arg("binwidth"), nb::arg("firstFrame"), nb::arg("finalFrame"));
}
