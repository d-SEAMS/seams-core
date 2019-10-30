#include <franzblau.hpp>

/********************************************/ /**
                                                *  Get all possible rings
                                                ***********************************************/
std::vector<std::vector<int>> primitive::countAllRings(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> neighHbondList, int maxDepth) {
  //
  primitive::Graph fullGraph;  // Contains all the information of the pointCloud
  std::vector<std::vector<int>> rings;  // List of all possible rings
  int nVertex;  // The number of vertices or points in the graph (in yCloud)

  // ------------------------------
  // Init
  // Initialize the graph object with all the information from the pointCloud
  // and neighbour list
  fullGraph = primitive::populateGraph(yCloud, neighHbondList);
  // ------------------------------
  // // Loop through every vertex
  // for (int iatom = 0; iatom < yCloud->nop; iatom++) {
  //   //
  //   // Recursive call
  // }  // loop through
  // ------------------------------

  return rings;
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
