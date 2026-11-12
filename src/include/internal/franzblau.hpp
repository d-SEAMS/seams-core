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
//! hydrogen-bonded neighbour list). Each ring starts at its lowest-indexed
//! member and runs toward the lesser of that member's ring neighbours, and
//! the rings are sorted lexicographically.
std::vector<std::vector<int>> ringNetwork(const std::vector<std::vector<int>> &nList,
                                           int maxDepth);

//! Rings on a CSR of cloud indices (no self entry, one edge each).
std::vector<std::vector<int>> ringNetwork(const std::vector<int> &offsets,
                                           const std::vector<int> &cols, int maxDepth);

//! The same ring set, built with a bounded ball at every vertex and a
//! primitivity test of every pair. ringNetwork leaves that table out when
//! the largest ring is at most 7, because a strict shortcut is then a walk
//! of at most two edges.
std::vector<std::vector<int>>
ringNetworkBalls(const std::vector<std::vector<int>> &nList, int maxDepth);

//! 0 when the last ringNetwork call answered every primitivity query from
//! the adjacency rows, 1 when it filled the sorted pair table, 2 when it
//! filled the occupancy bitset of the ball (rings of size 8 to 10).
[[nodiscard]] int lastRingSearchUsedBalls();

//! The rings ringNetwork returns whose lowest-indexed member is flagged in
//! @a sources (by index; empty flags every vertex), in the same order
std::vector<std::vector<int>> ringNetwork(const std::vector<std::vector<int>> &nList,
                                          int maxDepth,
                                          const std::vector<char> &sources);

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
 *  topology. A ring of at most maxDepth members can appear or vanish between
 *  frames only if it passes within (maxDepth/2 - 2)/2 bonds of a vertex
 *  whose neighbour row changed, which is on the ring itself up to maxDepth
 *  7. On a new frame the stored rings that meet that set are dropped and the
 *  rings through it are enumerated afresh. The output equals ringNetwork
 *  element by element, at a cost that scales with the number of changed
 *  rows rather than with system size.
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

  //! Vertices the last update enumerated rings through; the full vertex
  //! count on a whole-frame pass, zero when nothing changed
  [[nodiscard]] int lastRecomputedSources() const;

  //! Vertices whose bounded balls were rebuilt on the last update
  [[nodiscard]] int lastBallsRefreshed() const;

  //! Undirected bonds whose endpoints differ from the previous frame.
  //! Zero on the first frame and when the bond set is unchanged. The
  //! fraction below is this count divided by the molecule count, maintained
  //! while the rows are compared and readable in constant time afterwards.
  [[nodiscard]] int lastChangedBonds() const;

  //! Changed bonds per molecule on the last update.
  [[nodiscard]] double lastChangedBondFraction() const;

  //! Above this fraction of vertices with a changed row, update recomputes
  //! the whole frame in one pass instead. The default does not use a
  //! fraction. The update enumerates the sources in T and no others, so
  //! it takes a full pass only when T is the whole vertex set.
  //! One never falls back; zero always does. The result equals a rebuild
  //! at every fraction.
  void setFullPassFraction(double fraction);

  //! Above this many changed bonds per molecule, update takes the same
  //! whole-frame pass. The default is negative, which disables the test.
  //! The vertex fraction is then the only switch. Ignored when the
  //! vertex fraction is at least 1.
  void setBondCrossover(double changedBondsPerMolecule);

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

struct RingSearchSplit {
  double balls_ms = 0;
  double enumerate_ms = 0;
  double merge_ms = 0;
  double adj_ms = 0;
  double ball_count_ms = 0;
  double ball_fill_ms = 0;
};

//! Balls, source enumeration and the gather from the last ringNetwork call.
[[nodiscard]] RingSearchSplit lastRingSearchSplit();

} // namespace primitive

#endif // SEAMS_FRANZBLAU_H_
