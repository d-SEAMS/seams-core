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

#include <algorithm>
#include <vector>

#include <franzblau.hpp>

namespace {

/**
 * @brief Sorted (vertex, distance) pairs within a fixed radius of each vertex.
 * @details One bounded breadth-first sweep per vertex for the whole ring
 *  search, cleared by walking what each sweep touched. Every shortcut question
 *  the primitivity test asks is then a sorted lookup, in place of a fresh
 *  sweep per candidate ring.
 */
struct BoundedBalls {
  std::vector<std::vector<std::pair<int, int>>> ball;

  BoundedBalls(const std::vector<std::vector<int>> &adjacency, int radius)
      : ball(adjacency.size()) {
    const int nVertices = static_cast<int>(adjacency.size());
    std::vector<int> dist(nVertices, -1);
    std::vector<int> touched, frontier, next;
    for (int v = 0; v < nVertices; v++) {
      for (const int w : touched) {
        dist[w] = -1;
      }
      touched.clear();
      dist[v] = 0;
      touched.push_back(v);
      frontier.assign(1, v);
      for (int depth = 1; depth <= radius && !frontier.empty(); depth++) {
        next.clear();
        for (const int u : frontier) {
          for (const int w : adjacency[u]) {
            if (w >= 0 && w < nVertices && dist[w] == -1) {
              dist[w] = depth;
              touched.push_back(w);
              next.push_back(w);
            }
          }
        }
        frontier.swap(next);
      }
      auto &b = ball[v];
      b.reserve(touched.size());
      for (const int w : touched) {
        if (w != v) {
          b.emplace_back(w, dist[w]);
        }
      }
      std::sort(b.begin(), b.end());
    }
  }

  //! Graph distance from v to w when within the radius, else a large sentinel
  int dist(int v, int w) const {
    const auto &b = ball[v];
    auto it = std::lower_bound(b.begin(), b.end(), std::make_pair(w, 0));
    if (it != b.end() && it->first == w) {
      return it->second;
    }
    return 1 << 20;
  }
};

//! Bounded breadth-first levels from src, scratch reused across sources
void levelsFrom(const std::vector<std::vector<int>> &adjacency, int src,
                int maxLvl, std::vector<int> &lvl, std::vector<int> &touched,
                std::vector<int> &frontier, std::vector<int> &next) {
  for (const int w : touched) {
    lvl[w] = -1;
  }
  touched.clear();
  lvl[src] = 0;
  touched.push_back(src);
  frontier.assign(1, src);
  for (int depth = 1; depth <= maxLvl && !frontier.empty(); depth++) {
    next.clear();
    for (const int u : frontier) {
      for (const int w : adjacency[u]) {
        if (w >= 0 && w < static_cast<int>(lvl.size()) && lvl[w] == -1) {
          lvl[w] = depth;
          touched.push_back(w);
          next.push_back(w);
        }
      }
    }
    frontier.swap(next);
  }
}

//! All shortest paths src -> target under the level field, excluding src
void allShortestPaths(const std::vector<std::vector<int>> &adjacency,
                      const std::vector<int> &lvl, int src, int target,
                      std::vector<int> cur,
                      std::vector<std::vector<int>> &out) {
  if (target == src) {
    std::reverse(cur.begin(), cur.end());
    out.push_back(std::move(cur));
    return;
  }
  for (const int w : adjacency[target]) {
    if (w >= 0 && w < static_cast<int>(lvl.size()) && lvl[w] >= 0 &&
        lvl[w] == lvl[target] - 1) {
      std::vector<int> nxt = cur;
      nxt.push_back(target);
      allShortestPaths(adjacency, lvl, src, w, std::move(nxt), out);
    }
  }
}

//! Shortest-path criterion over every pair of members, by ball lookup
bool ringIsPrimitive(const BoundedBalls &balls, const std::vector<int> &ring) {
  const int k = static_cast<int>(ring.size());
  for (int i = 0; i < k; i++) {
    for (int j = i + 2; j < k; j++) {
      const int sep = j - i;
      const int ringDist = std::min(sep, k - sep);
      if (ringDist <= 1) {
        continue;
      }
      if (balls.dist(ring[i], ring[j]) < ringDist) {
        return false;
      }
    }
  }
  return true;
}

//! Append v if absent; ring sizes are tiny, so a scan beats a set
inline bool pushUniqueMember(std::vector<int> &ring, int v) {
  for (const int x : ring) {
    if (x == v) {
      return false;
    }
  }
  ring.push_back(v);
  return true;
}


/**
 * @brief Reusable scratch for the per-source ring enumeration.
 */
struct RingScratch {
  std::vector<int> lvl;
  std::vector<int> touched, frontier, next;
  std::vector<std::vector<int>> paths, pathsQ;
};

/**
 * @details Enumerates every primitive ring whose lowest-indexed member is
 *  @a src, appending to @a out. Shared by the whole-frame construction and
 *  the incremental updater, which re-runs it only for sources a
 *  frame-to-frame change can reach.
 */
void enumerateFromSource(const std::vector<std::vector<int>> &adjacency,
                         const BoundedBalls &balls, int src, int maxDepth,
                         int maxLvl, RingScratch &scr,
                         std::vector<std::vector<int>> &out) {
  const int nVertices = static_cast<int>(adjacency.size());
  if (static_cast<int>(scr.lvl.size()) != nVertices) {
    scr.lvl.assign(nVertices, -1);
    scr.touched.clear();
  }
  levelsFrom(adjacency, src, maxLvl, scr.lvl, scr.touched, scr.frontier,
             scr.next);
  // Directing: only the lowest-indexed member of a ring enumerates it, so
  // vertices below the source leave the level field entirely
  for (int v = 0; v < src; v++) {
    if (scr.lvl[v] >= 0) {
      scr.lvl[v] = -1;
    }
  }
  for (const int p : scr.touched) {
    if (scr.lvl[p] < 1 || scr.lvl[p] > maxLvl) {
      continue;
    }
    scr.paths.clear();
    allShortestPaths(adjacency, scr.lvl, src, p, {}, scr.paths);
    // Even rings: two vertex-disjoint shortest paths to an antipodal vertex
    if (2 * scr.lvl[p] >= 3 && 2 * scr.lvl[p] <= maxDepth) {
      for (size_t a = 0; a < scr.paths.size(); a++) {
        for (size_t b = a + 1; b < scr.paths.size(); b++) {
          std::vector<int> ring{src};
          bool ok = true;
          for (const int v : scr.paths[a]) {
            ok = ok && pushUniqueMember(ring, v);
          }
          for (int i = static_cast<int>(scr.paths[b].size()) - 2; i >= 0 && ok;
               i--) {
            ok = pushUniqueMember(ring, scr.paths[b][i]);
          }
          if (ok && ringIsPrimitive(balls, ring)) {
            out.push_back(std::move(ring));
          }
        }
      }
    }
    // Odd rings: an antipodal edge between two vertices at the same level
    if (2 * scr.lvl[p] + 1 <= maxDepth) {
      for (const int q : adjacency[p]) {
        if (q <= p || q >= nVertices || scr.lvl[q] != scr.lvl[p]) {
          continue;
        }
        scr.pathsQ.clear();
        allShortestPaths(adjacency, scr.lvl, src, q, {}, scr.pathsQ);
        for (const auto &pa : scr.paths) {
          for (const auto &pb : scr.pathsQ) {
            std::vector<int> ring{src};
            bool ok = true;
            for (const int v : pa) {
              ok = ok && pushUniqueMember(ring, v);
            }
            for (int i = static_cast<int>(pb.size()) - 1; i >= 0 && ok; i--) {
              ok = pushUniqueMember(ring, pb[i]);
            }
            if (ok && ringIsPrimitive(balls, ring)) {
              out.push_back(std::move(ring));
            }
          }
        }
      }
    }
  }
}

} // namespace

/**
 * @details Returns the primitive (shortest-path) rings, by index, given a
 *  neighbour list (also by index) and the largest ring size sought.
 *
 *  The construction follows Yuan and Cormack (Comput. Mater. Sci. 24, 343,
 *  2002): every member of a primitive ring sits at a graph distance from any
 *  other member equal to their separation around the ring, so each arc of the
 *  ring is a shortest path from any member taken as a source. A ring of even
 *  size @f$2L@f$ is therefore a source vertex plus two vertex-disjoint
 *  shortest paths to a vertex at level @f$L@f$; a ring of odd size
 *  @f$2L+1@f$ closes instead over an edge joining two vertices both at level
 *  @f$L@f$. Candidates are built only from shortest paths, and each ring is
 *  enumerated exactly once by restricting the search to its lowest-indexed
 *  member, so the output carries no duplicates by construction.
 *
 *  Candidates are then kept or dropped by the same criterion Franzblau's
 *  filter applied: no pair of members may be joined through the graph by
 *  fewer edges than around the ring. Each such question is answered by a
 *  lookup into per-vertex bounded neighbourhoods computed once for the whole
 *  search.
 *
 *  primitive::countAllRingsFromIndex and primitive::removeNonSPrings remain
 *  available and produce the same ring set by the generate-then-filter route.
 * @param[in] nList Row-ordered neighbour list by index (and NOT the atom ID)
 * @param[in] maxDepth The largest ring size to search for. Rings larger than
 *  maxDepth are not generated.
 * @return A vector of vectors of the rings; each ring contains the atom
 *  indices of the ring members.
 */
std::vector<std::vector<int>>
primitive::ringNetwork(const std::vector<std::vector<int>> &nList, int maxDepth) {
  std::vector<std::vector<int>> rings;
  if (maxDepth < 3 || nList.empty()) {
    return rings;
  }

  // Adjacency by index; the first element of each nList row is the vertex
  // itself, so it is skipped
  const int nVertices = static_cast<int>(nList.size());
  std::vector<std::vector<int>> adjacency(nVertices);
  for (int i = 0; i < nVertices; i++) {
    for (size_t j = 1; j < nList[i].size(); j++) {
      adjacency[i].push_back(nList[i][j]);
    }
  }

  const int maxLvl = maxDepth / 2;
  const int radius = std::max(maxLvl - 1, 1);
  const BoundedBalls balls(adjacency, radius);

  RingScratch scratch;
  for (int src = 0; src < nVertices; src++) {
    enumerateFromSource(adjacency, balls, src, maxDepth, maxLvl, scratch,
                        rings);
  }

  return rings;
}

/**
 *  @details Get all possible rings (only atom indices, not IDs). The input
 *   neighbour list is in terms of indices. All possible rings (including
 * non-SP) rings are generated using a recursive backtracking algorithm. This
 * function calls the following functions internally:
 *   - primitive::populateGraphFromIndices (for initializing the Graph object)
 *   - primitive::findRings (for getting all rings, using backtracking)
 *   - primitive::restoreEdgesFromIndices (restores back-up of the edges since
 * some may have been removed)
 *  @param[in] neighHbondList Row-ordered neighbour list by atom index (not ID).
 *  @param[in] maxDepth The maximum depth upto which rings will be searched.
 *   This means that rings larger than maxDepth will not be generated.
 *  @return The Graph object for the current frame.
 */
primitive::Graph
primitive::countAllRingsFromIndex(const std::vector<std::vector<int>> &neighHbondList,
                                  int maxDepth) {
  //
  primitive::Graph fullGraph; // Graph object
  std::vector<int>
      visited; // Contains the indices (NOT IDs) visited for a particular node
  int depth;   // Current depth

  // ------------------------------
  // Init
  // Initialize the graph object with all the information from the neighbour
  // list (of indices only)
  fullGraph = primitive::populateGraphFromIndices(neighHbondList);
  // ------------------------------
  // Loop through every vertex
  for (int iatom = 0; iatom < neighHbondList.size(); iatom++) {
    visited.clear();
    depth = 0;
    // Call the function for getting rings
    primitive::findRings(fullGraph, iatom, visited, maxDepth, depth);
  } // loop through every vertex
  // ------------------------------
  // Restore back-up of the edges (some may have been removed)
  primitive::restoreEdgesFromIndices(fullGraph, neighHbondList);
  // ------------------------------

  return fullGraph;
}

/**
 * @details All possible rings are searched for in this function, which
 * recursively calls itself. The rings are 'grown' from the root node (which is
 * the first vertex) using the backtracking algorithm. When it is first called
 * (before the root node has been assigned), root is a dummy value (which is
 * equal to -1, a value that can never be legitimate).
 *  @param[in, out] fullGraph Graph object containing the vertices (and the
 * neighbour lists). Vertices may be 'removed' from the Graph.
 *  @param[in] v The current vertex being visited or checked. It is added to the
 * list of all vertices visited.
 *  @param[in] visited A vector containing a list of the vertices visited for
 * book-keeping. If the current visited vector fulfills the condition for being
 * a ring, it is added to the rings vector of vector in the Graph.
 *  @param[in] maxDepth The maximum depth upto which rings will be searched.
 * This means that rings larger than maxDepth will not be generated.
 *  @param[in] depth The current depth. When this function is called for the
 * first time from primitive::countAllRingsFromIndex, the depth is initialized
 * to 0. When the depth is greater than or equal to maxDepth, the function
 * exits.
 *  @param[in] root The first vertex, from which the current visited vector
 * (candidate ring) is being grown. This is initialized to a dummy value of -1,
 * when it is called from primitive::countAllRingsFromIndex.
 */
int primitive::findRings(Graph &fullGraph, int v, std::vector<int> &visited,
                         int maxDepth, int depth, int root) {
  //
  int nnumNeighbours; // Number of nearest neighbours for iNode
  int n;              // Node of a neighbour (index)
  // -------------
  // For the first call:
  if (root == -1) {
    root = v;                          // Init the root as the current node
    fullGraph.pts[v].inGraph = false; // Remove the root from the graph
  }                                    // first call
  //
  // Checks
  if (depth >= maxDepth) {
    return 1; // false?
  }           // searched until the maximum length specified
  //
  // Add the current node to the visited vector
  visited.push_back(v);
  // -------------
  // Depth-first search
  // 1. Pick a root node (above)
  // 2. Search all the neighbouring nodes
  // 3. Start a recursion call, from the second nearest neighbours
  // 3(i) If the neighbour equals the root (selected in 1), ring has been found!
  // 3(ii) Or else search for all the unvisited neighbours
  // 4. Remove the root node from the graph
  // 5. Update the vector of vectors containing all rings
  // -------------
  // Start the search
  depth += 1; // Go to the next layer
  nnumNeighbours = fullGraph.pts[v].neighListIndex.size();
  // Loop through the neighbours of iNode
  for (int j = 0; j < nnumNeighbours; j++) {
    n = fullGraph.pts[v].neighListIndex[j]; // Neighbour index
    // Has a ring been found?!
    if (depth > 2 && n == root) {
      // Add the visited vector to the rings vector of vector
      fullGraph.rings.push_back(visited);
    } // A ring has been found!
    // Otherwise search all the neighbours which have not been searched already
    else if (fullGraph.pts[n].inGraph) {
      fullGraph.pts[n].inGraph = false; // Set to false now
      // Recursive call
      primitive::findRings(fullGraph, n, visited, maxDepth, depth, root);
      fullGraph.pts[n].inGraph = true; // Put n back in the graph
    } // Search the other unsearched neighbours
  }   // end of search through all neighbours
  // -------------
  // When the depth is 2, remove just the root from the neighbours of iNode
  if (depth == 2) {
    // Search for root in the neighbour list of v
    for (int j = 0; j < fullGraph.pts[v].neighListIndex.size(); j++) {
      if (root == fullGraph.pts[v].neighListIndex[j]) {
        fullGraph.pts[v].neighListIndex.erase(
            fullGraph.pts[v].neighListIndex.begin() + j);
      } // end of erase
    }   // end of search
  }     // remove root not edges

  //
  visited.pop_back();
  //
  return 0;
}

/**
 * @details Fills a Graph object with information from the PointCloud and the
 *  neighbour list. The indices in the neighbour list in the Vertex object are
 *  NOT the atom IDs (they are the atom indices according to the input
 *  PointCloud). The input neighbour list is by atom ID.
 * @param[in] yCloud The input PointCloud.
 * @param[in] neighHbondList The row-ordered neighbour list, containing atom
 *  IDs, and not the atom indices.
 * @return The Graph object for the current frame.
 */
primitive::Graph primitive::populateGraphFromNListID(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &neighHbondList) {
  //
  primitive::Graph fullGraph; // Contains all the information of the pointCloud
  primitive::Vertex iVertex;  // The vertex corresponding to a particular point
  int nnumNeighbours;         // Number of nearest neighbours for iatom
  std::vector<int> iNeigh;    // Neighbours of the current iatom
  int jatomID;                // Atom ID of the nearest neighbour
  int jatomIndex;             // Atom index of the nearest neighbour
  // ------------------------------
  // Loop through every point in yCloud
  for (int iatom = 0; iatom < yCloud.nop; iatom++) {
    // -----
    // Update the neighbours of iatom (only the indices!)
    iNeigh.clear(); // init
    nnumNeighbours = neighHbondList[iatom].size() -
                     1; // The first element is atomID of iatom
    for (int j = 1; j <= nnumNeighbours; j++) {
      jatomID = neighHbondList[iatom][j]; // Atom ID
      // Get the atom index for the vector nearest neighbour list
      auto it = yCloud.idIndexMap.find(jatomID);
      if (it != yCloud.idIndexMap.end()) {
        jatomIndex = it->second;
      } // found jatomIndex
      else {
        std::cerr << "Panic induced!\n";
        return fullGraph;
      } // jatomIndex not found
      // Update iNeigh
      iNeigh.push_back(jatomIndex);
    } // end of loop through nearest neighbours
    // -----
    // Update iVertex
    iVertex.atomIndex = iatom;
    iVertex.neighListIndex = iNeigh;
    // Add to the Graph object
    fullGraph.pts.push_back(iVertex);
  } // end of loop through every iatom
  // ------------------------------

  return fullGraph;
}

/**
 * @details Fills a Graph object with information from the PointCloud and the
 *  neighbour list. The indices in the neighbour list in the Vertex object are
 *  NOT the atom IDs (they are the atom indices according to the input
 *  PointCloud). The input neighbour list is by index NOT atom IDs. Otherwise,
 *  this function does the same thing as primitive::populateGraphFromNListID.
 * The only difference is that this function takes the neighbour list BY INDEX.
 * @param[in] nList The row-ordered neighbour list, containing atom
 *  indices (according to the input PointCloud).
 * @return The Graph object for the current frame.
 */
primitive::Graph
primitive::populateGraphFromIndices(const std::vector<std::vector<int>> &nList) {
  //
  primitive::Graph fullGraph; // Contains all the information of the pointCloud
  primitive::Vertex iVertex;  // The vertex corresponding to a particular point
  int nnumNeighbours;         // Number of nearest neighbours for iatom
  int iatom, jatom;           // Atom index being saved
  // ------------------------------
  // Loop through every point in nList
  for (int i = 0; i < nList.size(); i++) {
    iatom = nList[i][0]; // Atom index of i
    // neighListIndex is simply the i^th row of nList
    //
    // Update iVertex
    iVertex.atomIndex = iatom;
    iVertex.neighListIndex = nList[i];
    // Add to the Graph object
    fullGraph.pts.push_back(iVertex);
  } // end of loop through iatom

  return fullGraph;
}

/**
 * @details Re-fills the neighbour lists of a Graph object from a row-ordered
 *  neighbour list (which is BY INDEX not IDs). Some vertices may have been
 *  removed while rings were generated using the backtracking algorithm
 *  (primitive::findRings). Also, the indices in the neighbour list in the
 *  Vertex object are not the atom IDs.
 * @param[in] fullGraph The Graph object for the current frame. The neighbour
 *  lists of component Vertex objects may have been depleted.
 * @param[in] nList The row-ordered neighbour list, containing atom
 *  indices (according to the input PointCloud).
 * @return The Graph object for the current frame.
 */
void
primitive::restoreEdgesFromIndices(Graph &fullGraph,
                                   const std::vector<std::vector<int>> &nList) {
  //
  // ------------------------------
  // Loop through every point in nList
  for (int i = 0; i < nList.size(); i++) {
    // neighListIndex is simply the i^th row of nList
    // Update the neighListIndex list in the graph object
    fullGraph.pts[i].neighListIndex = nList[i];
  } // end of loop through iatom
}

/**
 * @details Removes non-SP rings (according to the Franzblau shortest path
 *  criterion) from the rings vector of vectors member n the Graph object. The
 *  rings vector of vectors contains indices and NOT atom IDs. This function
 *  calls primitive::shortestPath internally to calculate the shortest path.
 * @param[in] fullGraph The Graph object for the current frame. This also
 *  contains the rings vector of vectors, which has all possible rings (possibly
 *  inclding non-SP rings).
 * @return The Graph object for the current frame.
 */
void primitive::removeNonSPrings(primitive::Graph &fullGraph) {
  //
  const int nVertices = fullGraph.pts.size(); // Number of vertices in the graph
  const int nRings = fullGraph.rings.size();  // Number of rings
  std::vector<bool> ringsToRemove; // Vector containing the logical values for
                                   // removal of the current ring index
  std::vector<std::vector<int>>
      emptyTempRings; // Empty vector of vectors to swap
  std::vector<std::vector<int>>
      primitiveRings; // Vector of vectors of rings after removing non SP rings
  // -------------------
  // Make sure all the vertices are in the graph before removing non-SP rings
  for (int iVer = 0; iVer < nVertices; iVer++) {
    fullGraph.pts[iVer].inGraph = true;
  } // end of loop through every vertex
  // -------------------
  // The criterion compares, for each pair of non-adjacent ring members, the
  // separation around the ring against the separation through the graph. Only
  // the latter needs searching, it is bounded by half the longest ring, and it
  // does not depend on which ring the pair was drawn from.
  //
  // Collect every such question first and group them by the vertex they start
  // from. One breadth-first sweep per distinct starting vertex then answers all
  // of its questions at once, in place of one depth-first search per pair per
  // ring. The sweep reaches only the handful of vertices within the bound, so
  // its scratch is cleared by walking what it touched rather than the whole
  // graph.
  ringsToRemove.assign(nRings, false);

  struct Query {
    int target;   // The other ring member
    int ringHops; // Separation around the ring, in edges
    int ring;     // Which ring asked
  };
  std::vector<std::vector<Query>> queries(nVertices);

  int maxHops = 0;
  for (int iRing = 0; iRing < nRings; iRing++) {
    const std::vector<int> &currentRing = fullGraph.rings[iRing];
    const int ringSize = currentRing.size();
    maxHops = std::max(maxHops, ringSize / 2);
    for (int jVer = 0; jVer < ringSize; jVer++) {
      // connect with all other, skip j-j (distance=0) and j-(j+1) (nearest
      // neighbors)
      for (int kVer = jVer + 2; kVer < ringSize; kVer++) {
        const int currentV = currentRing[jVer];
        const int currentN = currentRing[kVer];
        if (currentV < 0 || currentV >= nVertices || currentN < 0 ||
            currentN >= nVertices) {
          continue;
        }
        const int d_jk = std::abs(jVer - kVer);
        queries[currentV].push_back(
            {currentN, std::min(d_jk, std::abs(d_jk - ringSize)), iRing});
      } // end of loop through k^th vertex
    }   // end of loop through j^th vertex
  }     // end of loop through rings

  // Answer each vertex's questions with one bounded sweep
  std::vector<int> distance(nVertices, -1);
  std::vector<int> touched;
  std::vector<int> frontier;
  std::vector<int> next;

  for (int v = 0; v < nVertices; v++) {
    if (queries[v].empty()) {
      continue;
    }

    touched.clear();
    frontier.assign(1, v);
    distance[v] = 0;
    touched.push_back(v);
    for (int depth = 1; depth <= maxHops && !frontier.empty(); depth++) {
      next.clear();
      for (const int u : frontier) {
        for (const int w : fullGraph.pts[u].neighListIndex) {
          if (w >= 0 && w < nVertices && distance[w] == -1) {
            distance[w] = depth;
            touched.push_back(w);
            next.push_back(w);
          }
        }
      }
      frontier.swap(next);
    }

    for (const Query &q : queries[v]) {
      // Separation through the graph, in edges; -1 when beyond the bound
      const int graphHops = distance[q.target];
      if (graphHops >= 0 && graphHops < q.ringHops) {
        ringsToRemove[q.ring] = true;
      } // Decide whether to keep or remove the ring
    }

    // Clear only what this sweep reached
    for (const int w : touched) {
      distance[w] = -1;
    }
  } // end of loop through vertices
  // -------------------
  // Remove all the rings whose indices are given in the ringsToRemove vector
  for (int i = 0; i < ringsToRemove.size(); i++) {
    if (!ringsToRemove[i]) {
      primitiveRings.push_back(fullGraph.rings[i]);
    } // updates new copy
  }   // end of loop through ringsToRemove
  // -------------------
  // Update the graph rings with the primitiveRings
  emptyTempRings.swap(fullGraph.rings);
  fullGraph.rings = std::move(primitiveRings);
  // -------------------
}

/**
 * @details Calculates the shortest path for a particular ring. This function
 * uses recursion.
 * @param[in] fullGraph The Graph object for the current frame.
 * @param[in] v The current vertex being checked.
 * @param[in] goal The first element of the candidate ring being checked (the
 *  root node).
 * @param[in] path The path or length of the visited points (This basically
 *  contains all the indices in the visited vector, excluding the current
 * vertex).
 * @param[in] visited This vector contains the indices of the vertices visited
 *  or checked (for book-keeping).
 * @param[in] maxDepth The maximum depth or maximum length of the rings.
 * @param[in] depth The current depth. When this function is called from
 *  primitive::removeNonSPrings, the depth is initialized as 0.
 * @return The Graph object for the current frame.
 */
int primitive::shortestPath(Graph &fullGraph, int v, int goal,
                            std::vector<int> &path, std::vector<int> &visited,
                            int maxDepth, int depth) {
  int len_path = 0;   // Length of the path
  int nnumNeighbours; // Number of neighbours for a particular v
  int n;              // Index of the nearest neighbour of v
  // Start the search for the shortest path
  if (depth < maxDepth) {
    depth += 1; // One layer below
    visited.push_back(
        v); // Add the current vertex to the path (visited points)
    //
    if (v == goal) {
      len_path = path.size(); // Path of the length of visited points
      // If the current path is shorter OR this is the first path found
      if (depth < len_path || len_path == 0) {
        path = visited;
        maxDepth = depth;
      } // Current path is the shortest
    }   // Goal reached
    // Recursive calls to function
    else {
      nnumNeighbours = fullGraph.pts[v].neighListIndex.size();
      // Search all the neighbours
      for (int j = 0; j < nnumNeighbours; j++) {
        n = fullGraph.pts[v].neighListIndex[j]; // Index of nearest neighbour
        // If n has not already been searched:
        if (fullGraph.pts[n].inGraph == true) {
          fullGraph.pts[n].inGraph = false; // Set to false
          primitive::shortestPath(fullGraph, n, goal, path, visited, maxDepth,
                                  depth);   // Recursive call
          fullGraph.pts[n].inGraph = true; // Put back in the graph
        } // If n is in the graph, call recursively
      }   // End of loop over all neighbours
    }     // Goal not reached
    //
    // Pop the vector
    visited.pop_back();
  } // for depth less than maxDepth
  return 0;
}

/**
 * @details Function for clearing Graph if it is already
 *  filled. This should be called before every frame is read in.
 * @param[out] currentGraph The cleared Graph
 */
primitive::Graph primitive::clearGraph(Graph &currentGraph) {
  //
  std::vector<primitive::Vertex> tempPts;
  std::vector<std::vector<int>> tempRings;
  tempPts.swap(currentGraph.pts);
  tempRings.swap(currentGraph.rings);
  return currentGraph;
}

/**
 * @details Implementation state of primitive::RingUpdater.
 *
 *  The locality bound that makes the update exact: every member of a ring of
 *  size at most maxDepth lies within maxLvl = maxDepth/2 hops of the ring's
 *  lowest-indexed member, and the primitivity of the ring is decided by paths
 *  of fewer than maxLvl edges between members. A source further than
 *  2*maxLvl + 1 hops from every endpoint of a changed edge, measured in the
 *  union of the old and the new graph, therefore encloses its rings and every
 *  path that could decide them entirely in unchanged territory: its ring set
 *  cannot differ between frames. Distances in the union graph never exceed
 *  those in either frame's graph, so measuring there errs only toward
 *  recomputing more.
 */
struct primitive::RingUpdater::Impl {
  int maxDepth = 0;
  int maxLvl = 0;
  int ballRadius = 1;
  std::vector<std::vector<int>> adjacency;              // previous frame
  std::vector<std::vector<std::vector<int>>> bySource;  // rings, by source
  std::vector<std::vector<int>> flattened;              // last returned set
  int recomputed = 0;
};

primitive::RingUpdater::RingUpdater(int maxDepth)
    : impl_(new Impl) {
  impl_->maxDepth = maxDepth;
  impl_->maxLvl = maxDepth / 2;
  impl_->ballRadius = std::max(impl_->maxLvl - 1, 1);
}

primitive::RingUpdater::~RingUpdater() = default;
primitive::RingUpdater::RingUpdater(RingUpdater &&) noexcept = default;
primitive::RingUpdater &
primitive::RingUpdater::operator=(RingUpdater &&) noexcept = default;

int primitive::RingUpdater::lastRecomputedSources() const {
  return impl_->recomputed;
}

const std::vector<std::vector<int>> &
primitive::RingUpdater::update(const std::vector<std::vector<int>> &nList) {
  Impl &st = *impl_;
  const int nVertices = static_cast<int>(nList.size());

  // Adjacency rows sorted, so frame-to-frame comparison is a vector compare
  std::vector<std::vector<int>> newAdj(nVertices);
  for (int i = 0; i < nVertices; i++) {
    for (size_t j = 1; j < nList[i].size(); j++) {
      newAdj[i].push_back(nList[i][j]);
    }
    std::sort(newAdj[i].begin(), newAdj[i].end());
  }

  const bool fullPass =
      st.adjacency.size() != static_cast<size_t>(nVertices) ||
      st.maxDepth < 3;

  if (fullPass) {
    st.adjacency = std::move(newAdj);
    st.bySource.assign(nVertices, {});
    // Balls rebuild in linear time each frame; the incrementality that pays
    // is confined to the enumeration
    RingScratch scratch;
    const BoundedBalls balls(st.adjacency, st.ballRadius);
    st.recomputed = nVertices;
    for (int src = 0; src < nVertices; src++) {
      enumerateFromSource(st.adjacency, balls, src, st.maxDepth, st.maxLvl,
                          scratch, st.bySource[src]);
    }
  } else {
    // Vertices whose adjacency row changed: the endpoints of every added or
    // removed edge
    std::vector<int> changed;
    for (int v = 0; v < nVertices; v++) {
      if (st.adjacency[v] != newAdj[v]) {
        changed.push_back(v);
      }
    }

    if (!changed.empty()) {
      // Union graph reaches everything either frame can reach
      std::vector<std::vector<int>> unionAdj(nVertices);
      for (int v = 0; v < nVertices; v++) {
        unionAdj[v] = st.adjacency[v];
        unionAdj[v].insert(unionAdj[v].end(), newAdj[v].begin(),
                           newAdj[v].end());
        std::sort(unionAdj[v].begin(), unionAdj[v].end());
        unionAdj[v].erase(
            std::unique(unionAdj[v].begin(), unionAdj[v].end()),
            unionAdj[v].end());
      }

      const int reach = 2 * st.maxLvl - 1;
      std::vector<int> dist(nVertices, -1);
      std::vector<int> frontier, next;
      for (const int v : changed) {
        dist[v] = 0;
        frontier.push_back(v);
      }
      std::vector<int> affected = changed;
      for (int depth = 1; depth <= reach && !frontier.empty(); depth++) {
        next.clear();
        for (const int u : frontier) {
          for (const int w : unionAdj[u]) {
            if (w >= 0 && w < nVertices && dist[w] == -1) {
              dist[w] = depth;
              affected.push_back(w);
              next.push_back(w);
            }
          }
        }
        frontier.swap(next);
      }

      st.adjacency = std::move(newAdj);

      const BoundedBalls balls(st.adjacency, st.ballRadius);
      RingScratch scratch;
      st.recomputed = 0;
      for (const int src : affected) {
        st.bySource[src].clear();
        enumerateFromSource(st.adjacency, balls, src, st.maxDepth, st.maxLvl,
                            scratch, st.bySource[src]);
        st.recomputed++;
      }
    } else {
      st.recomputed = 0;
    }
  }

  st.flattened.clear();
  for (const auto &group : st.bySource) {
    st.flattened.insert(st.flattened.end(), group.begin(), group.end());
  }
  return st.flattened;
}
