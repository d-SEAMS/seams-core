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
  int atomID;                       // This may not be the same as the index
  double x, y, z;                   // The coordinates of the vertex
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
Graph populateGraph(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                    std::vector<std::vector<int>> neighHbondList);

// Creates a vector of vectors of all possible rings
std::vector<std::vector<int>> countAllRings(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> neighHbondList, int maxDepth);

}  // namespace primitive

#endif  // __FRANZBLAU_H_
