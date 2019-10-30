#include <franzblau.hpp>

/********************************************/ /**
                                                *  Get all possible rings (only
                                                *atom indices, not IDs)
                                                ***********************************************/
primitive::Graph primitive::countAllRings(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> neighHbondList, int maxDepth) {
  //
  primitive::Graph fullGraph;  // Contains all the information of the pointCloud
  std::vector<std::vector<int>>
      rings;  // List of all possible rings (update at the end) (SAVE ATOM IDs)
  std::vector<int>
      visited;  // Contains the indices (NOT IDs) visited for a particular node
  int depth;    // Current depth
  // Variables required for restoring the neighbour list
  std::vector<int> iNeigh;  // Neighbour list of indices NOT ATOM IDs
  int nnumNeighbours;       // Number of nearest neighbours for iatom
  int jatomID, jatomIndex;  // For the neighbour

  // ------------------------------
  // Init
  // Initialize the graph object with all the information from the pointCloud
  // and neighbour list
  fullGraph = primitive::populateGraph(yCloud, neighHbondList);
  // ------------------------------
  // Loop through every vertex
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    // init
    visited.clear();
    depth = 0;  // zero for the first call
    // Recursive call
    primitive::findRings(&fullGraph, iatom, &visited, maxDepth, &depth);
  }  // loop through
  // ------------------------------
  // Restore back-up of the edges (some may have been removed)
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
        std::cerr << "Panic induced! Something is very wrong with your map.\n";
      }  // jatomIndex not found
      // Update iNeigh
      iNeigh.push_back(jatomIndex);
    }  // end of loop through nearest neighbours
    // -----
    // Update fullGraph
    fullGraph.pts[iatom].neighListIndex = iNeigh;
  }  // end of loop through every iatom
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
int primitive::findRings(Graph *fullGraph, int iNode, std::vector<int> *visited,
                         int maxDepth, int *depth, int root) {
  //
  int nnumNeighbours;  // Number of nearest neighbours for iNode
  int n;               // Node of a neighbour (index)
  // -------------
  // For the first call:
  if (root == -1) {
    root = iNode;                           // Init the root as the current node
    fullGraph->pts[iNode].inGraph = false;  // Remove the root from the graph
  }                                         // first call
  //
  // Checks
  if (*depth >= maxDepth) {
    return 1;  // false?
  }            // searched until the maximum length specified
  //
  // Add the current node to the visited vector
  visited->push_back(iNode);
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
  *depth += 1;  // Go to the next layer
  nnumNeighbours = fullGraph->pts[iNode].neighListIndex.size();
  // Loop through the neighbours of iNode
  for (int j = 0; j < nnumNeighbours; j++) {
    n = fullGraph->pts[iNode].neighListIndex[j];  // Neighbour index
    // Has a ring been found?!
    if (*depth > 2 && n == root) {
      // Add the visited vector to the rings vector of vector
      fullGraph->rings.push_back((*visited));
      continue;
    }  // A ring has been found!
    // Otherwise search all the neighbours which have not been searched already
    if (fullGraph->pts[n].inGraph) {
      fullGraph->pts[n].inGraph = false;  // Set to false now
      // Recursive call
      primitive::findRings(fullGraph, n, visited, maxDepth, depth, root);
      fullGraph->pts[n].inGraph = true;  // Put n back in the graph
    }  // Search the other unsearched neighbours
  }    // end of search through all neighbours
  // -------------
  // When the depth is 2, remove just the root from the neighbours of iNode
  if (*depth == 2) {
    // Search for root in the neighbour list of iNode
    auto it = std::find(fullGraph->pts[iNode].neighListIndex.begin(),
                        fullGraph->pts[iNode].neighListIndex.end(), root);
    // If found, remove root
    if (it != fullGraph->pts[iNode].neighListIndex.end()) {
      fullGraph->pts[iNode].neighListIndex.erase(it);
    }  // search and remove in neighbour list
  }    // remove root not edges

  //
  (*visited).pop_back();
  //
  return 0;
}

/********************************************/ /**
                                                *  Fills a Graph object with
                                                *information from pointCloud and
                                                *the neighbour list. The indices
                                                * in the neighbour list in the
                                                *Vertex object are NOT the atom
                                                *IDs.
                                                ***********************************************/
primitive::Graph primitive::populateGraph(
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
    iVertex.atomID = yCloud->pts[iatom].atomID;
    iVertex.x = yCloud->pts[iatom].x;  // x coordinate
    iVertex.y = yCloud->pts[iatom].y;  // y coordinate
    iVertex.z = yCloud->pts[iatom].z;  // z coordinate
    iVertex.neighListIndex = iNeigh;
    // Add to the Graph object
    fullGraph.pts.push_back(iVertex);
  }  // end of loop through every iatom
  // ------------------------------

  return fullGraph;
}
