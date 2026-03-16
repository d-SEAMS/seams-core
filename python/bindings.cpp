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

    nb::class_<molSys::Point<double>>(m, "PointDouble",
        "Per-particle data: coordinates, type, molecule ID, ice classification.")
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

    nb::enum_<molSys::bond_type>(m, "BondType",
        "Bond classification: staggered, eclipsed, or out_of_range.")
        .value("staggered", molSys::bond_type::staggered)
        .value("eclipsed", molSys::bond_type::eclipsed)
        .value("out_of_range", molSys::bond_type::out_of_range);

    nb::enum_<molSys::atom_state_type>(m, "AtomStateType",
        "Per-atom ice phase classification from CHILL/CHILL+/q6.")
        .value("cubic", molSys::atom_state_type::cubic)
        .value("hexagonal", molSys::atom_state_type::hexagonal)
        .value("water", molSys::atom_state_type::water)
        .value("interfacial", molSys::atom_state_type::interfacial)
        .value("clathrate", molSys::atom_state_type::clathrate)
        .value("interClathrate", molSys::atom_state_type::interClathrate)
        .value("unclassified", molSys::atom_state_type::unclassified)
        .value("reCubic", molSys::atom_state_type::reCubic)
        .value("reHex", molSys::atom_state_type::reHex);

    nb::class_<molSys::Result>(m, "Result",
        "Bond correlation result: classifier (bond type) and c_value.")
        .def(nb::init<>())
        .def_rw("classifier", &molSys::Result::classifier)
        .def_rw("c_value", &molSys::Result::c_value)
        .def("__repr__", [](const molSys::Result &self_C) {
            std::uintptr_t ptr_val = std::uintptr_t(&self_C);
            return std::format("<Result mem_loc:{:x}>", static_cast<uint64_t>(ptr_val));
        });

    nb::class_<molSys::PointCloud<molSys::Point<double>, double>>(m, "PointCloudDouble",
        "Collection of points for a single frame, with box dimensions.")
        .def(nb::init<>())
        .def_rw("pts", &molSys::PointCloud<molSys::Point<double>, double>::pts)
        .def_rw("currentFrame", &molSys::PointCloud<molSys::Point<double>, double>::currentFrame)
        .def_rw("nop", &molSys::PointCloud<molSys::Point<double>, double>::nop)
        .def_rw("box", &molSys::PointCloud<molSys::Point<double>, double>::box)
        .def_rw("boxLow", &molSys::PointCloud<molSys::Point<double>, double>::boxLow)
        .def_rw("idIndexMap", &molSys::PointCloud<molSys::Point<double>, double>::idIndexMap);

    // I/O (lambdas hide the yCloud* in/out parameter, creating it internally)
    m.def("readXYZ", &sinp::readXYZ,
          "Read atom coordinates from an XYZ file.",
          nb::arg("filename"));
    m.def("readLammpsTrjreduced",
          [](std::string filename, int targetFrame, int typeI, bool isSlice,
             std::array<double, 3> coordLow, std::array<double, 3> coordHigh) {
            molSys::PointCloud<molSys::Point<double>, double> yCloud;
            return sinp::readLammpsTrjreduced(filename, targetFrame, yCloud,
                                               typeI, isSlice, coordLow, coordHigh);
          },
          "Read a LAMMPS trajectory frame, keeping only atoms of the given type.",
          nb::arg("filename"), nb::arg("targetFrame"), nb::arg("typeI"),
          nb::arg("isSlice"), nb::arg("coordLow"), nb::arg("coordHigh"));
    m.def("readLammpsTrjO",
          [](std::string filename, int targetFrame, int typeO, bool isSlice,
             std::array<double, 3> coordLow, std::array<double, 3> coordHigh) {
            molSys::PointCloud<molSys::Point<double>, double> yCloud;
            return sinp::readLammpsTrjO(filename, targetFrame, yCloud,
                                         typeO, isSlice, coordLow, coordHigh);
          },
          "Read a LAMMPS trajectory frame, keeping only oxygen atoms.",
          nb::arg("filename"), nb::arg("targetFrame"), nb::arg("typeO"),
          nb::arg("isSlice"), nb::arg("coordLow"), nb::arg("coordHigh"));
    m.def("readLammpsTrj",
          [](std::string filename, int targetFrame, bool isSlice,
             std::array<double, 3> coordLow, std::array<double, 3> coordHigh) {
            molSys::PointCloud<molSys::Point<double>, double> yCloud;
            return sinp::readLammpsTrj(filename, targetFrame, yCloud,
                                        isSlice, coordLow, coordHigh);
          },
          "Read a LAMMPS trajectory frame with all atom types.",
          nb::arg("filename"), nb::arg("targetFrame"),
          nb::arg("isSlice"), nb::arg("coordLow"), nb::arg("coordHigh"));
    m.def("readBonds", &sinp::readBonds,
          "Read bond connectivity from a formatted bond file.",
          nb::arg("filename"));
    m.def("atomInSlice", &sinp::atomInSlice,
          "Check whether a point (x, y, z) lies within a volume slice.",
          nb::arg("x"), nb::arg("y"), nb::arg("z"),
          nb::arg("coordLow"), nb::arg("coordHigh"));

    // Neighbours
    m.def("clearNeighbourList", &nneigh::clearNeighbourList,
          "Free memory for a neighbour list.",
          nb::arg("nList"));
    m.def("getNewNeighbourListByIndex", &nneigh::getNewNeighbourListByIndex,
          "Build a neighbour list by index using a distance cutoff.",
          nb::arg("yCloud"), nb::arg("cutoff"));
    m.def("halfNeighList", &nneigh::halfNeighList,
          "Build a half neighbour list (each pair stored once) for one atom type.",
          nb::arg("yCloud"), nb::arg("rcutoff"), nb::arg("typeI"));
    m.def("neighbourListByIndex", &nneigh::neighbourListByIndex,
          "Convert an atom-ID neighbour list to an index-based neighbour list.",
          nb::arg("yCloud"), nb::arg("nList"));
    m.def("neighList", &nneigh::neighList,
          "Build a full neighbour list for two atom types within a cutoff.",
          nb::arg("yCloud"), nb::arg("rcutoff"), nb::arg("typeI"), nb::arg("typeJ"));
    m.def("neighListO", &nneigh::neighListO,
          "Build a full neighbour list for a single atom type within a cutoff.",
          nb::arg("rcutoff"), nb::arg("yCloud"), nb::arg("typeI"));

    // Bonds
    m.def("createBondsFromCages", &bond::createBondsFromCages,
          "Create bond connectivity from rings and cage information.",
          nb::arg("rings"), nb::arg("cageList"), nb::arg("type"), nb::arg("nRings"));
    m.def("getHbondDistanceOH", &bond::getHbondDistanceOH,
          "Compute the O-H hydrogen bond distance between two atoms.",
          nb::arg("oCloud"), nb::arg("hCloud"), nb::arg("oAtomIndex"), nb::arg("hAtomIndex"));
    m.def("populateHbonds", &bond::populateHbonds,
          "Build the hydrogen bond network from a trajectory and neighbour list.",
          nb::arg("filename"), nb::arg("yCloud"), nb::arg("nList"),
          nb::arg("targetFrame"), nb::arg("Htype"));
    m.def("populateHbondsWithInputClouds", &bond::populateHbondsWithInputClouds,
          "Build hydrogen bonds using pre-loaded oxygen and hydrogen point clouds.",
          nb::arg("yCloud"), nb::arg("hCloud"), nb::arg("nList"));
    m.def("trimBonds", &bond::trimBonds,
          "Remove duplicate bonds from a bond list.",
          nb::arg("bonds"));

    // Ring analysis (Franzblau)
    m.def("clearGraph", &primitive::clearGraph,
          "Free memory for a graph object.",
          nb::arg("currentGraph"));
    m.def("countAllRingsFromIndex", &primitive::countAllRingsFromIndex,
          "Find all possible rings (including non-shortest-path) up to maxDepth.",
          nb::arg("neighHbondList"), nb::arg("maxDepth"));
    m.def("ringNetwork", &primitive::ringNetwork,
          "Find all primitive (shortest-path) rings up to maxDepth.",
          nb::arg("nList"), nb::arg("maxDepth"));
    m.def("populateGraphFromIndices", &primitive::populateGraphFromIndices,
          "Create a graph object from an index-based neighbour list.",
          nb::arg("nList"));
    m.def("populateGraphFromNListID", &primitive::populateGraphFromNListID,
          "Create a graph object from an atom-ID neighbour list and point cloud.",
          nb::arg("yCloud"), nb::arg("neighHbondList"));
    m.def("removeNonSPrings", &primitive::removeNonSPrings,
          "Remove non-shortest-path rings using the Franzblau criterion.",
          nb::arg("fullGraph"));
    m.def("restoreEdgesFromIndices", &primitive::restoreEdgesFromIndices,
          "Restore graph edges from an index-based neighbour list.",
          nb::arg("fullGraph"), nb::arg("nList"));

    // Ring classification
    m.def("assignPolygonType", &ring::assignPolygonType,
          "Assign atom types based on the ring size of n-membered rings.",
          nb::arg("rings"), nb::arg("atomTypes"), nb::arg("nRings"));
    m.def("assignPrismType", &ring::assignPrismType,
          "Assign atom types for atoms belonging to prism rings.",
          nb::arg("rings"), nb::arg("listPrism"), nb::arg("ringSize"),
          nb::arg("ringType"), nb::arg("atomTypes"), nb::arg("atomState"));
    m.def("clearRingList", &ring::clearRingList,
          "Free memory for a list of rings.",
          nb::arg("rings"));
    m.def("compareRings", &ring::compareRings,
          "Check whether two unordered rings contain the same elements.",
          nb::arg("ring1"), nb::arg("ring2"));
    m.def("commonElementsInThreeRings", &ring::commonElementsInThreeRings,
          "Check whether three rings share at least one common element.",
          nb::arg("ring1"), nb::arg("ring2"), nb::arg("ring3"));
    m.def("deformedPrismTypes", &ring::deformedPrismTypes,
          "Get atom type values for deformed prisms.",
          nb::arg("atomState"), nb::arg("atomTypes"), nb::arg("maxDepth"));
    m.def("discardExtraTetragonBlocks", &ring::discardExtraTetragonBlocks,
          "Discard duplicate 4-membered ring pairs that are parallel in one dimension.",
          nb::arg("basal1"), nb::arg("basal2"), nb::arg("yCloud"));
    m.def("findPrisms", &ring::findPrisms,
          "Identify which rings form prism blocks.",
          nb::arg("rings"), nb::arg("ringType"), nb::arg("nPerfectPrisms"),
          nb::arg("nImperfectPrisms"), nb::arg("nList"), nb::arg("rmsdPerAtom"),
          nb::arg("doShapeMatching"), nb::arg("yCloud"));
    m.def("findsCommonElements", &ring::findsCommonElements,
          "Return the common elements shared by two rings.",
          nb::arg("ring1"), nb::arg("ring2"));
    m.def("findTripletInRing", &ring::findTripletInRing,
          "Search for a triplet of atoms within a ring.",
          nb::arg("ring"), nb::arg("triplet"));
    m.def("getSingleRingSize", &ring::getSingleRingSize,
          "Extract rings of a specific size from a list of all rings.",
          nb::arg("rings"), nb::arg("ringSize"));
    m.def("hasCommonElements", &ring::hasCommonElements,
          "Check whether two rings share any common elements.",
          nb::arg("ring1"), nb::arg("ring2"));
    m.def("basalPrismConditions", &ring::basalPrismConditions,
          "Test whether two rings satisfy strict basal prism conditions.",
          nb::arg("nList"), nb::arg("basal1"), nb::arg("basal2"));
    m.def("relaxedPrismConditions", &ring::relaxedPrismConditions,
          "Test whether two rings satisfy relaxed prism conditions (at least one bond).",
          nb::arg("nList"), nb::arg("basal1"), nb::arg("basal2"));

    // Topology analysis
    m.def("polygonRingAnalysis", &ring::polygonRingAnalysis,
          "Classify rings in a quasi-2D monolayer system and write output.",
          nb::arg("path"), nb::arg("rings"), nb::arg("nList"), nb::arg("yCloud"),
          nb::arg("maxDepth"), nb::arg("sheetArea"), nb::arg("firstFrame"));
    m.def("bulkPolygonRingAnalysis", &ring::bulkPolygonRingAnalysis,
          "Classify rings in a bulk system and write output.",
          nb::arg("path"), nb::arg("rings"), nb::arg("nList"), nb::arg("yCloud"),
          nb::arg("maxDepth"), nb::arg("firstFrame"));
    m.def("prismAnalysis", &ring::prismAnalysis,
          "Run prism identification on rings up to maxDepth and write output.",
          nb::arg("path"), nb::arg("rings"), nb::arg("nList"), nb::arg("yCloud"),
          nb::arg("maxDepth"), nb::arg("atomID"), nb::arg("firstFrame"),
          nb::arg("currentFrame"), nb::arg("doShapeMatching"));
    m.def("rmAxialTranslations", &ring::rmAxialTranslations,
          "Remove axial translations from an ice nanotube for visualization.",
          nb::arg("yCloud"), nb::arg("atomID"), nb::arg("firstFrame"), nb::arg("currentFrame"));

    // TUM (Topological Unit Matching)
    m.def("atomsFromCages", &tum3::atomsFromCages,
          "Get atom indices belonging to cages in a given cluster.",
          nb::arg("rings"), nb::arg("cageList"), nb::arg("clusterCages"));
    m.def("averageRMSDatom", &tum3::averageRMSDatom,
          "Average the per-atom RMSD over the number of shared cages.",
          nb::arg("rmsdPerAtom"), nb::arg("noOfCommonAtoms"));
    m.def("buildRefDDC", &tum3::buildRefDDC,
          "Build a reference double-diamond cage from a template XYZ file.",
          nb::arg("fileName"));
    m.def("buildRefHC", &tum3::buildRefHC,
          "Build a reference hexagonal cage from a template XYZ file.",
          nb::arg("fileName"));
    m.def("clusterCages", &tum3::clusterCages,
          "Cluster cages using Stillinger's algorithm and write XYZ output.",
          nb::arg("yCloud"), nb::arg("path"), nb::arg("rings"), nb::arg("cageList"),
          nb::arg("numHC"), nb::arg("numDDC"));
    m.def("shapeMatchDDC", &tum3::shapeMatchDDC,
          "Shape-match a target double-diamond cage against a reference.",
          nb::arg("yCloud"), nb::arg("refPoints"), nb::arg("cageList"),
          nb::arg("cageIndex"), nb::arg("rings"), nb::arg("quat"), nb::arg("rmsd"));
    m.def("shapeMatchHC", &tum3::shapeMatchHC,
          "Shape-match a target hexagonal cage against a reference.",
          nb::arg("yCloud"), nb::arg("refPoints"), nb::arg("cageUnit"),
          nb::arg("rings"), nb::arg("nList"), nb::arg("quat"), nb::arg("rmsd"));
    m.def("topoBulkCriteria", &tum3::topoBulkCriteria,
          "Find HCs and DDCs in a bulk system using topological criteria.",
          nb::arg("path"), nb::arg("rings"), nb::arg("nList"), nb::arg("yCloud"),
          nb::arg("firstFrame"), nb::arg("numHC"), nb::arg("numDDC"), nb::arg("ringType"));
    m.def("topoUnitMatchingBulk", &tum3::topoUnitMatchingBulk,
          "Run full topological unit matching for bulk water.",
          nb::arg("path"), nb::arg("rings"), nb::arg("nList"), nb::arg("yCloud"),
          nb::arg("firstFrame"), nb::arg("printClusters"), nb::arg("onlyTetrahedral"));
    m.def("updateRMSDatom", &tum3::updateRMSDatom,
          "Update per-atom RMSD from a cage shape-matching result.",
          nb::arg("rings"), nb::arg("cageUnit"), nb::arg("rmsd"),
          nb::arg("rmsdPerAtom"), nb::arg("noOfCommonAtoms"), nb::arg("atomTypes"));

    // Selection
    m.def("getPointCloudOneAtomType", &gen::getPointCloudOneAtomType,
          "Extract a point cloud containing only atoms of a given type.",
          nb::arg("yCloud"), nb::arg("outCloud"), nb::arg("atomTypeI"),
          nb::arg("isSlice"), nb::arg("coordLow"), nb::arg("coordHigh"));
    m.def("atomsInSingleSlice", &gen::atomsInSingleSlice,
          "Mark atoms inside a rectangular volume slice.",
          nb::arg("yCloud"), nb::arg("clearPreviousSliceSelection"),
          nb::arg("coordLow"), nb::arg("coordHigh"));
    m.def("moleculesInSingleSlice", &gen::moleculesInSingleSlice,
          "Mark whole molecules as in-slice if any atom is inside the region.",
          nb::arg("yCloud"), nb::arg("clearPreviousSliceSelection"),
          nb::arg("coordLow"), nb::arg("coordHigh"));
    m.def("setAtomsWithSameMolID", &gen::setAtomsWithSameMolID,
          "Set the inSlice flag for all atoms sharing a given molecule ID.",
          nb::arg("yCloud"), nb::arg("molIDAtomIDmap"), nb::arg("molID"), nb::arg("inSliceValue"));
    m.def("getEdgeMoleculesInRings", &ring::getEdgeMoleculesInRings,
          "Select edge molecules in rings that straddle the slice boundary.",
          nb::arg("rings"), nb::arg("oCloud"), nb::arg("yCloud"),
          nb::arg("identicalCloud"), nb::arg("coordLow"), nb::arg("coordHigh"));
    m.def("printSliceGetEdgeMoleculesInRings", &ring::printSliceGetEdgeMoleculesInRings,
          "Select edge molecules in rings and write slice output files.",
          nb::arg("path"), nb::arg("rings"), nb::arg("oCloud"), nb::arg("yCloud"),
          nb::arg("coordLow"), nb::arg("coordHigh"), nb::arg("identicalCloud"));

    // CHILL/CHILL+ classification
    m.def("getCorrelPlus", &chill::getCorrelPlus,
          "Compute CHILL+ bond-order correlations and classify bond types.",
          nb::arg("yCloud"), nb::arg("nList"), nb::arg("isSlice"));
    m.def("getIceTypePlus", &chill::getIceTypePlus,
          "Classify each atom's ice type using CHILL+ and write to file.",
          nb::arg("yCloud"), nb::arg("nList"), nb::arg("path"),
          nb::arg("firstFrame"), nb::arg("isSlice"), nb::arg("outputFileName"));
    m.def("getCorrel", &chill::getCorrel,
          "Compute CHILL bond-order correlations and classify bond types.",
          nb::arg("yCloud"), nb::arg("nList"), nb::arg("isSlice"));
    m.def("getIceType", &chill::getIceType,
          "Classify each atom's ice type using CHILL and write to file.",
          nb::arg("yCloud"), nb::arg("nList"), nb::arg("path"),
          nb::arg("firstFrame"), nb::arg("isSlice"), nb::arg("outputFileName"));
    m.def("getq6", &chill::getq6,
          "Compute the q6 bond order parameter for all atoms.",
          nb::arg("yCloud"), nb::arg("nList"), nb::arg("isSlice"));
    m.def("reclassifyWater", &chill::reclassifyWater,
          "Reclassify water molecules using averaged q6 and q3 parameters.",
          nb::arg("yCloud"), nb::arg("q6"));
    m.def("printIceType", &chill::printIceType,
          "Print the ice type classification for the current frame.",
          nb::arg("yCloud"), nb::arg("path"), nb::arg("firstFrame"),
          nb::arg("isSlice"), nb::arg("outputFileName"));

    // Output
    m.def("writeDump", &sout::writeDump,
          "Write a LAMMPS dump file for the current point cloud.",
          nb::arg("yCloud"), nb::arg("path"), nb::arg("outFile"));

    // Clustering
    m.def("clusterAnalysis", &clump::clusterAnalysis,
          "Cluster ice-like particles and return the largest ice cluster.",
          nb::arg("path"), nb::arg("iceCloud"), nb::arg("yCloud"), nb::arg("nList"),
          nb::arg("iceNeighbourList"), nb::arg("cutoff"), nb::arg("firstFrame"),
          nb::arg("bopAnalysis"));
    m.def("recenterClusterCloud", &clump::recenterClusterCloud,
          "Recenter a cluster point cloud for visualization.",
          nb::arg("iceCloud"), nb::arg("nList"));

    // RDF
    m.def("rdf2Danalysis_AA", &rdf2::rdf2Danalysis_AA,
          "Compute the 2D radial distribution function for identical atom types.",
          nb::arg("path"), nb::arg("rdfValues"), nb::arg("yCloud"),
          nb::arg("cutoff"), nb::arg("binwidth"), nb::arg("firstFrame"), nb::arg("finalFrame"));
}
