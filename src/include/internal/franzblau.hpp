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

#ifndef SEAMS_FRANZBLAU_H_
#define SEAMS_FRANZBLAU_H_

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cmath>
#include <memory>
#include <sstream>
#include <string>
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
 * Franzblau shortest path criterion (primitive::removeNonSPrings), answered
 * by a hop-bounded breadth-first sweep per vertex.
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
std::vector<std::vector<int>> ringNetwork(const std::vector<std::vector<int>> &nList,
                                          int maxDepth);

//! Creates a graph object and fills it with the information from a neighbour
//! list and pointCloud created before. NOTE: the neighbourListIndex contains
//! the indices and NOT the atom IDs as in the neighbour list
Graph populateGraphFromNListID(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &neighHbondList);

//! Creates a graph object and fills it with the information from a neighbour
//! list of INDICES NOT ATOM IDs created before. NOTE: the neighbourListIndex
//! contains the indices and NOT the atom IDs as in the neighbour list
Graph populateGraphFromIndices(const std::vector<std::vector<int>> &nList);

//! Re-fills the neighbour lists of a graph object from a neighbour
//! list of INDICES NOT ATOM IDs created before. NOTE: the neighbourListIndex
//! contains the indices and NOT the atom IDs as in the neighbour list
void restoreEdgesFromIndices(Graph &fullGraph,
                             const std::vector<std::vector<int>> &nList);

//! Creates a vector of vectors of all possible rings
Graph countAllRingsFromIndex(const std::vector<std::vector<int>> &neighHbondList,
                             int maxDepth);

//! Removes the non-SP rings, using the Franzblau shortest path criterion
void removeNonSPrings(Graph &fullGraph);

//! Main function that searches for all rings
[[nodiscard]] int findRings(Graph &fullGraph, int v, std::vector<int> &visited, int maxDepth,
              int depth, int root = -1);

//! Calculates the shortest path
[[nodiscard]] int shortestPath(Graph &fullGraph, int v, int goal, std::vector<int> &path,
                 std::vector<int> &visited, int maxDepth, int depth = 1);

//! Function for clearing vectors in Graph after multiple usage
Graph clearGraph(Graph &currentGraph);


/** @class primitive::RingUpdater
 * @brief Exact frame-to-frame maintenance of the primitive ring network.
 * @details Consecutive trajectory frames share almost all of their bond
 *  topology. Rings are stored partitioned by their lowest-indexed member, and
 *  on a new frame only sources within a proven locality radius of a changed
 *  edge are re-enumerated: a ring's members and every path that can decide
 *  its primitivity lie within a bounded neighbourhood of its source, so a
 *  source far enough from every change keeps its rings unchanged. The output
 *  equals a full recomputation exactly, at a cost that scales with the edge
 *  churn between frames rather than with system size.
 */
class RingUpdater {
public:
  explicit RingUpdater(int maxDepth);
  ~RingUpdater();
  RingUpdater(RingUpdater &&) noexcept;
  RingUpdater &operator=(RingUpdater &&) noexcept;
  RingUpdater(const RingUpdater &) = delete;
  RingUpdater &operator=(const RingUpdater &) = delete;

  //! Primitive rings for this frame's neighbour list (row-ordered, by index,
  //! first element of each row the vertex itself)
  const std::vector<std::vector<int>> &
  update(const std::vector<std::vector<int>> &nList);

  //! Sources re-enumerated by the last update; the full vertex count on the
  //! first frame, zero when nothing changed
  [[nodiscard]] int lastRecomputedSources() const;

  //! Vertices whose bounded balls were rebuilt on the last update
  [[nodiscard]] int lastBallsRefreshed() const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace primitive

#endif // SEAMS_FRANZBLAU_H_
