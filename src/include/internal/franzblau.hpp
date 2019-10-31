#ifndef __FRANZBLAU_H_
#define __FRANZBLAU_H_

#include <math.h>
#include <sys/stat.h>
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <cage.hpp>
#include <mol_sys.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>

namespace primitive {

// A vertex is a collection of elements required for graph traversal
struct Vertex {
  int atomIndex;                    // This is the index according to pointCloud
  std::vector<int> neighListIndex;  // contains the INDICES (not the atomIDs) of
                                    // the neighbouring vertices
  bool inGraph =
      true;  // True by default. Setting it to false removes it from the graph
  bool visited =
      false;  // Flag determining whether the vertex has been visited or not
};

// A graph is a per-frame object, containing all the vertices and a vector of
// vectors of the rings found
struct Graph {
  std::vector<Vertex> pts;  // Collection of vertices. The index of each should
                            // be the same as that in pointCloud
  std::vector<std::vector<int>>
      rings;  // List of all the rings (of every size) found
};

// Creates a graph object and fills it with the information from a neighbour
// list and pointCloud created before. NOTE: the neighbourListIndex contains the
// indices and NOT the atom IDs as in the neighbour list
Graph populateGraphFromNListID(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> neighHbondList);

// Creates a graph object and fills it with the information from a neighbour
// list of INDICES NOT ATOM IDs created before. NOTE: the neighbourListIndex
// contains the indices and NOT the atom IDs as in the neighbour list
Graph populateGraphFromIndices(std::vector<std::vector<int>> nList);

// Re-fills the neighbour lists of a graph object from a neighbour
// list of INDICES NOT ATOM IDs created before. NOTE: the neighbourListIndex
// contains the indices and NOT the atom IDs as in the neighbour list
Graph restoreEdgesFromIndices(Graph *fullGraph,
                              std::vector<std::vector<int>> nList);

// Creates a vector of vectors of all possible rings
Graph countAllRingsFromIndex(std::vector<std::vector<int>> neighHbondList,
                             int maxDepth);

// Creates a vector of vectors of all possible rings
Graph removeNonSPrings(Graph *fullGraph);

// Main function that searches for all rings
int findRings(Graph *fullGraph, int iNode, std::vector<int> *visited,
              int maxDepth, int depth, int root = -1);

// Calculates the shortest path
int shortestPath(Graph *fullGraph, int v, int goal, std::vector<int> *path,
                 std::vector<int> *visited, int maxDepth, int depth = 1);

//// Function for clearing vectors in Graph after multiple usage
Graph clearGraph(Graph *currentGraph);

}  // namespace primitive

#endif  // __FRANZBLAU_H_
