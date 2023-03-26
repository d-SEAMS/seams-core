//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#ifndef __FRANZBLAU_H_
#define __FRANZBLAU_H_

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <math.h>
#include <memory>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

#include <cage.hpp>
#include <mol_sys.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>

/** @file franzblau.hpp
 *    @brief File for generating shortest-path rings according to the Franzblau
 *   algorithm.
 */

/**
 *  @addtogroup primitive
 *  @{
 */

/** @brief Functions for generating primitive rings.
 *         This namespace contains struct definitions and functions that are
 * used for generating primitive (shortest-path) rings (directed cyclic graphs).
 *
 * The Vertex object is a collection of elements for each point, required for
 * graph traversal. The Graph object is an object for the whole frame,
 * containing the information of all vertices, and a row-ordered vector of
 * vector of the rings generated.
 *
 * The <a
 * href="https://journals.aps.org/prb/pdf/10.1103/PhysRevB.44.4925">Franzblau
 * shortest-path criterion</a> has been used. The SP (shortest-path) criterion
 * is midway between the least restrictive and most restrictive criteria in the
 * hierarchy.
 *
 * The following is the procedure for finding primitive rings:
 *
 * 1. All possible rings (including non-SP) rings are found, in the
 * primitive::countAllRingsFromIndex function, using the <a
 * href="https://www.math.upenn.edu/~wilf/website/CombinatorialAlgorithms.pdf">backtracking
 * algorithm</a>. This is a recursive algorithm.
 *
 * 2. The non-SP rings are then removed from the list of all rings, using the
 * Franzblau shortest path criterion (primitive::removeNonSPrings). This also
 * uses recursion.
 *
 *   ### Changelog ###
 *
 *  - Amrita Goswami [amrita16thaug646@gmail.com]; date modified: Nov 14, 2019
 *  - Rohit Goswami [rog32@hi.is]; date modified: Mar 20, 2021
 */

namespace primitive {

/** @struct Vertex
 * @brief This is a collection of elements, for each point, required for graph
 * traversal.
 *
 * Contains specifically the members:
 * - @b atomIndex : This is the index according to the PointCloud.
 * - @b neighListIndex : A vector of indices (not IDs) of the neighboring
 * vertices.
 * - @b inGraph : Bool qualifier, which is true by default. Setting it to
 * false removes the vertex from the graph.
 */
struct Vertex {
  int atomIndex;                   //! This is the index according to pointCloud
  std::vector<int> neighListIndex; //! Contains the INDICES (not the atomIDs)
                                   //! of the neighbouring vertices
  bool inGraph =
      true; //! True by default. Setting it to false removes it from the graph
};

/*! @struct Graph
 * @brief This is a per-frame object, containing all the vertices for the
 * particular frame, along with the vector of rings generated.
 *
 * Contains specifically the members:
 * - @b pts : Collection of vertices. The index of each should be the index
 * according to the PointCloud.
 * - @b rings : A row-ordered vector of vectors for the rings generated,
 * containing the indices (not IDs) of each member of the rings.
 */
struct Graph {
  std::vector<Vertex> pts; //! Collection of vertices. The index of each should
                           //! be the same as that in pointCloud
  std::vector<std::vector<int>>
      rings; //! List of all the rings (of every size) found
};

//! Returns a vector of vectors containing the rings (of all sizes), by atom
//! index, given the neighbour list also by index (preferably the
//! hydrogen-bonded neighbour list). Internally uses the Graph and Vertex
//! objects.
std::vector<std::vector<int>> ringNetwork(std::vector<std::vector<int>> nList,
                                          int maxDepth);

//! Creates a graph object and fills it with the information from a neighbour
//! list and pointCloud created before. NOTE: the neighbourListIndex contains
//! the indices and NOT the atom IDs as in the neighbour list
Graph populateGraphFromNListID(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> neighHbondList);

//! Creates a graph object and fills it with the information from a neighbour
//! list of INDICES NOT ATOM IDs created before. NOTE: the neighbourListIndex
//! contains the indices and NOT the atom IDs as in the neighbour list
Graph populateGraphFromIndices(std::vector<std::vector<int>> nList);

//! Re-fills the neighbour lists of a graph object from a neighbour
//! list of INDICES NOT ATOM IDs created before. NOTE: the neighbourListIndex
//! contains the indices and NOT the atom IDs as in the neighbour list
Graph restoreEdgesFromIndices(Graph *fullGraph,
                              std::vector<std::vector<int>> nList);

//! Creates a vector of vectors of all possible rings
Graph countAllRingsFromIndex(std::vector<std::vector<int>> neighHbondList,
                             int maxDepth);

//! Removes the non-SP rings, using the Franzblau shortest path criterion
Graph removeNonSPrings(Graph *fullGraph);

//! Main function that searches for all rings
int findRings(Graph *fullGraph, int v, std::vector<int> *visited, int maxDepth,
              int depth, int root = -1);

//! Calculates the shortest path
int shortestPath(Graph *fullGraph, int v, int goal, std::vector<int> *path,
                 std::vector<int> *visited, int maxDepth, int depth = 1);

//! Function for clearing vectors in Graph after multiple usage
Graph clearGraph(Graph *currentGraph);

} // namespace primitive

#endif // __FRANZBLAU_H_
