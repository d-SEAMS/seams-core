#include <franzblau.hpp>

/********************************************/ /**
                                                *  The vector of vector of
                                                *rings, by index, is returned,
                                                *given a neighbour list (also by
                                                *index) and the maximum depth
                                                *upto which rings will be
                                                *searched, using the Frnazblau
                                                *algorithm for shortest paths.
                                                *This function is registered in
                                                *Lua and exposed to the user.
                                                ***********************************************/
std::vector<std::vector<int>> primitive::ringNetwork(
    std::vector<std::vector<int>> nList, int maxDepth) {
  //
  primitive::Graph fullGraph;  // Graph object, contains the connectivity
                               // information from the neighbourlist

  // Find all possible rings, using backtracking. This may contain non-primitive
  // rings as well.
  fullGraph = primitive::countAllRingsFromIndex(nList, maxDepth);

  // Remove all non-SP rings using the Franzblau algorithm.
  fullGraph = primitive::removeNonSPrings(&fullGraph);

  // The rings vector of vectors inside the fullGraph graph object is the ring
  // network information we want
  return fullGraph.rings;
}

/********************************************/ /**
                                                *  Get all possible rings (only
                                                *atom indices, not IDs). The
                                                *neighbour list is in terms of
                                                *indices
                                                ***********************************************/
primitive::Graph primitive::countAllRingsFromIndex(
    std::vector<std::vector<int>> neighHbondList, int maxDepth) {
  //
  primitive::Graph fullGraph;  // Graph object
  std::vector<int>
      visited;  // Contains the indices (NOT IDs) visited for a particular node
  int depth;    // Current depth

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
    primitive::findRings(&fullGraph, iatom, &visited, maxDepth, depth);
  }  // loop through every vertex
  // ------------------------------
  // Restore back-up of the edges (some may have been removed)
  fullGraph = primitive::restoreEdgesFromIndices(&fullGraph, neighHbondList);
  // ------------------------------

  return fullGraph;
}

/********************************************/ /**
                                                *  All possible rings are
                                                *searched for here in this
                                                *function, which recursively
                                                *calls itself. When it is first
                                                *called, root is a dummy value
                                                *(-1)
                                                ***********************************************/
int primitive::findRings(Graph *fullGraph, int v, std::vector<int> *visited,
                         int maxDepth, int depth, int root) {
  //
  int nnumNeighbours;  // Number of nearest neighbours for iNode
  int n;               // Node of a neighbour (index)
  // -------------
  // For the first call:
  if (root == -1) {
    root = v;                           // Init the root as the current node
    fullGraph->pts[v].inGraph = false;  // Remove the root from the graph
  }                                     // first call
  //
  // Checks
  if (depth >= maxDepth) {
    return 1;  // false?
  }            // searched until the maximum length specified
  //
  // Add the current node to the visited vector
  visited->push_back(v);
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
  depth += 1;  // Go to the next layer
  nnumNeighbours = fullGraph->pts[v].neighListIndex.size();
  // Loop through the neighbours of iNode
  for (int j = 0; j < nnumNeighbours; j++) {
    n = fullGraph->pts[v].neighListIndex[j];  // Neighbour index
    // Has a ring been found?!
    if (depth > 2 && n == root) {
      // Add the visited vector to the rings vector of vector
      fullGraph->rings.push_back(*visited);
    }  // A ring has been found!
    // Otherwise search all the neighbours which have not been searched already
    else if (fullGraph->pts[n].inGraph) {
      fullGraph->pts[n].inGraph = false;  // Set to false now
      // Recursive call
      primitive::findRings(fullGraph, n, visited, maxDepth, depth, root);
      fullGraph->pts[n].inGraph = true;  // Put n back in the graph
    }  // Search the other unsearched neighbours
  }    // end of search through all neighbours
  // -------------
  // When the depth is 2, remove just the root from the neighbours of iNode
  if (depth == 2) {
    // Search for root in the neighbour list of v
    for (int j = 0; j < fullGraph->pts[v].neighListIndex.size(); j++) {
      if (root == fullGraph->pts[v].neighListIndex[j]) {
        fullGraph->pts[v].neighListIndex.erase(
            fullGraph->pts[v].neighListIndex.begin() + j);
      }  // end of erase
    }    // end of search
    // // Search for root in the neighbour list of v
    // auto it = std::find(fullGraph->pts[v].neighListIndex.begin(),
    //                     fullGraph->pts[v].neighListIndex.end(), root);
    // // If found, remove root
    // if (it != fullGraph->pts[v].neighListIndex.end()) {
    //   fullGraph->pts[v].neighListIndex.erase(it);
    // }  // search and remove in neighbour list
  }  // remove root not edges

  //
  (*visited).pop_back();
  // (*visited).erase(visited->end() - 1);
  //
  // return 0;
}

/********************************************/ /**
                                                *  Fills a Graph object with
                                                *information from pointCloud and
                                                *the neighbour list. The indices
                                                * in the neighbour list in the
                                                *Vertex object are NOT the atom
                                                *IDs.
                                                ***********************************************/
primitive::Graph primitive::populateGraphFromNListID(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> neighHbondList) {
  //
  primitive::Graph fullGraph;  // Contains all the information of the pointCloud
  primitive::Vertex iVertex;   // The vertex corresponding to a particular point
  int nnumNeighbours;          // Number of nearest neighbours for iatom
  std::vector<int> iNeigh;     // Neighbours of the current iatom
  int jatomID;                 // Atom ID of the nearest neighbour
  int jatomIndex;              // Atom index of the nearest neighbour
  // ------------------------------
  // Loop through every point in yCloud
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    // -----
    // Update the neighbours of iatom (only the indices!)
    iNeigh.clear();  // init
    nnumNeighbours = neighHbondList[iatom].size() -
                     1;  // The first element is atomID of iatom
    for (int j = 1; j <= nnumNeighbours; j++) {
      jatomID = neighHbondList[iatom][j];  // Atom ID
      // Get the atom index for the vector nearest neighbour list
      auto it = yCloud->idIndexMap.find(jatomID);
      if (it != yCloud->idIndexMap.end()) {
        jatomIndex = it->second;
      }  // found jatomIndex
      else {
        std::cerr << "Panic induced!\n";
        return fullGraph;
      }  // jatomIndex not found
      // Update iNeigh
      iNeigh.push_back(jatomIndex);
    }  // end of loop through nearest neighbours
    // -----
    // Update iVertex
    iVertex.atomIndex = iatom;
    iVertex.neighListIndex = iNeigh;
    // Add to the Graph object
    fullGraph.pts.push_back(iVertex);
  }  // end of loop through every iatom
  // ------------------------------

  return fullGraph;
}

/********************************************/ /**
                                                *  Fills a Graph object with
                                                *information from pointCloud and
                                                *the neighbour list of INDICES
                                                *not atom IDs. The indices in
                                                *the neighbour list in the
                                                *Vertex object are NOT the atom
                                                *IDs.
                                                ***********************************************/
primitive::Graph primitive::populateGraphFromIndices(
    std::vector<std::vector<int>> nList) {
  //
  primitive::Graph fullGraph;  // Contains all the information of the pointCloud
  primitive::Vertex iVertex;   // The vertex corresponding to a particular point
  int nnumNeighbours;          // Number of nearest neighbours for iatom
  int iatom, jatom;            // Atom index being saved
  // ------------------------------
  // Loop through every point in nList
  for (int i = 0; i < nList.size(); i++) {
    iatom = nList[i][0];  // Atom index of i
    // neighListIndex is simply the i^th row of nList
    //
    // Update iVertex
    iVertex.atomIndex = iatom;
    iVertex.neighListIndex = nList[i];
    // Add to the Graph object
    fullGraph.pts.push_back(iVertex);
  }  // end of loop through iatom

  return fullGraph;
}

/********************************************/ /**
                                                *  Re-fills the neighbour lists
                                                *of a graph object from a
                                                *neighbour the neighbour list of
                                                *INDICES not atom IDs. The
                                                *indices in the neighbour list
                                                *in the Vertex object are NOT
                                                *the atom IDs.
                                                ***********************************************/
primitive::Graph primitive::restoreEdgesFromIndices(
    Graph *fullGraph, std::vector<std::vector<int>> nList) {
  //
  // ------------------------------
  // Loop through every point in nList
  for (int i = 0; i < nList.size(); i++) {
    // neighListIndex is simply the i^th row of nList
    // Update the neighListIndex list in the graph object
    fullGraph->pts[i].neighListIndex = nList[i];
  }  // end of loop through iatom

  return *fullGraph;
}

/********************************************/ /**
                                                *  Removes non-SP rings
                                                *(Franzblau algorithm) from the
                                                *rings vector of vectors member
                                                *in the Graph object. The rings
                                                *vector still only contains
                                                *indices and NOT atom IDs.
                                                ***********************************************/
primitive::Graph primitive::removeNonSPrings(primitive::Graph *fullGraph) {
  //
  int nVertices = fullGraph->pts.size();  // Number of vertices in the graph
  int nRings = fullGraph->rings.size();   // Number of rings
  std::vector<bool> ringsToRemove;  // Vector containing the logical values for
                                    // removal of the current ring index
  std::vector<int> currentRing;     // Current ring being evaluated
  int ringSize;                     // Length of the current ring
  bool removeRing;  // Logical for removing the current ring (true) or not
  std::vector<std::vector<int>>
      emptyTempRings;  // Empty vector of vectors to swap
  std::vector<std::vector<int>>
      primitiveRings;  // Vector of vectors of rings after removing non SP rings
  int currentV;        // Current vertex
  int currentN;        // Current neighbour
  int dist_r;          // Distance over ring
  int dist_g;          // Distance over the entire graph
  int d_jk;            // Absolute difference between j and k
  std::vector<int> path;     // Vector containing a path
  std::vector<int> visited;  // Vector containing the visited points
  // -------------------
  // Make sure all the vertices are in the graph before removing non-SP rings
  for (int iVer = 0; iVer < nVertices; iVer++) {
    fullGraph->pts[iVer].inGraph = true;
  }  // end of loop through every vertex
  // -------------------
  // Loop through every ring
  for (int iRing = 0; iRing < nRings; iRing++) {
    currentRing = fullGraph->rings[iRing];  // Current ring
    ringSize = currentRing.size();          // Length of the current ring
    removeRing = false;                     // init
    // Loop through every j^th vertex
    for (int jVer = 0; jVer < ringSize; jVer++) {
      // connect with all other, skip j-j (distance=0) and j-(j+1) (nearest
      // neighbors)
      for (int kVer = jVer + 2; kVer < ringSize; kVer++) {
        // If not remove
        if (!removeRing) {
          currentV = currentRing[jVer];  // current 'vertex'
          currentN = currentRing[kVer];  // Current 'neighbour'
          d_jk = std::abs(jVer - kVer);
          dist_r = std::min(d_jk, std::abs(d_jk - ringSize)) +
                   1;    // Distance over the ring
          path.clear();  // init
          visited.clear();
          // Call shortest path function
          primitive::shortestPath(fullGraph, currentV, currentN, &path,
                                  &visited, dist_r, 0);
          dist_g = path.size();  // Length of the path over the graph
          if (dist_g < dist_r) {
            removeRing = true;
          }  // Decide whether to keep or remove the ring
        }    // If ring is not to be removed
      }      // end of loop through k^th vertex
    }        // end of loop through j^th vertex
    // Update bool value for removal of currentRing
    ringsToRemove.push_back(removeRing);
  }  // end of loop through rings
  // -------------------
  // Remove all the rings whose indices are given in the ringsToRemove vector
  for (int i = 0; i < ringsToRemove.size(); i++) {
    if (!ringsToRemove[i]) {
      primitiveRings.push_back(fullGraph->rings[i]);
    }  // updates new copy
  }    // end of loop through ringsToRemove
  // -------------------
  // Update the graph rings with the primitiveRings
  emptyTempRings.swap(fullGraph->rings);
  fullGraph->rings = primitiveRings;
  // -------------------
  return *fullGraph;
}

/********************************************/ /**
                                                *  Calculates the shortest path
                                                *for a particular ring
                                                ***********************************************/
int primitive::shortestPath(Graph *fullGraph, int v, int goal,
                            std::vector<int> *path, std::vector<int> *visited,
                            int maxDepth, int depth) {
  int len_path = 0;    // Length of the path
  int nnumNeighbours;  // Number of neighbours for a particular v
  int n;               // Index of the nearest neighbour of v
  // Start the search for the shortest path
  if (depth < maxDepth) {
    depth += 1;  // One layer below
    (*visited).push_back(
        v);  // Add the current vertex to the path (visited points)
    //
    if (v == goal) {
      len_path = (*path).size();  // Path of the length of visited points
      // If the current path is shorter OR this is the first path found
      if (depth < len_path || len_path == 0) {
        (*path) = (*visited);
        maxDepth = depth;
      }  // Current path is the shortest
    }    // Goal reached
    // Recursive calls to function
    else {
      nnumNeighbours = fullGraph->pts[v].neighListIndex.size();
      // Search all the neighbours
      for (int j = 0; j < nnumNeighbours; j++) {
        n = fullGraph->pts[v].neighListIndex[j];  // Index of nearest neighbour
        // If n has not already been searched:
        if (fullGraph->pts[n].inGraph == true) {
          fullGraph->pts[n].inGraph = false;  // Set to false
          primitive::shortestPath(fullGraph, n, goal, path, visited, maxDepth,
                                  depth);    // Recursive call
          fullGraph->pts[n].inGraph = true;  // Put back in the graph
        }  // If n is in the graph, call recursively
      }    // End of loop over all neighbours
    }      // Goal not reached
    //
    // Pop the vector
    (*visited).pop_back();
  }  // for depth less than maxDepth
  // return 0;
}

/********************************************/ /**
 *  Function for clearing Graph if it is already
 filled. This should be called before every frame is read in.
 *  @param[out] currentGraph The cleared Graph
 ***********************************************/
primitive::Graph primitive::clearGraph(Graph *currentGraph) {
  //
  std::vector<primitive::Vertex> tempPts;
  std::vector<std::vector<int>> tempRings;
  tempPts.swap(currentGraph->pts);
  tempRings.swap(currentGraph->rings);
  return *currentGraph;
}
