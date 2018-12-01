#include <neighbors.hpp>

/********************************************/ /**
 *  Constructor
 ***********************************************/
neigh::treeKNN::treeKNN() {
  // Initialize things
  this->frame = new CMolecularSystem;
  // Volume Limits
  this->coordLow = this->coordHigh = {0, 0, 0};
}
/********************************************/ /**
 *  Destructor
 ***********************************************/
neigh::treeKNN::~treeKNN() { delete frame; }

/********************************************/ /**
 *  Prepare the point cloud for KNN
 ***********************************************/
void neigh::treeKNN::prepFrame(int nop, std::string fileName) {
  this->frame->initializeFrames(nop, fileName);
}

void neigh::treeKNN::populateCloud(int typeI) {
  int filteredParticles(0);
  int nop(frame->parameter->nop);
  int dummy(0);
  // Do filtering
  for (int t; t < nop; t++) {
    // Match type
    if (frame->molecules[t].type == typeI) {
      // Check limits
      if (this->isThere(t, frame)) {
        // Accomodate one more point
        cloud.pts.resize(cloud.pts.size() + 1);
        // Dump point
        double coordX = frame->molecules[t].get_posx();
        double coordY = frame->molecules[t].get_posy();
        double coordZ = frame->molecules[t].get_posz();
        cloud.pts[dummy].x = coordX;
        cloud.pts[dummy].y = coordY;
        cloud.pts[dummy].z = coordZ;
        // Update dummy
        dummy++;
      }
    }
  }
  cloud.box[0] = frame->parameter->boxx;
  cloud.box[1] = frame->parameter->boxy;
  cloud.box[2] = frame->parameter->boxz;
}

// Implements knnSearch to get the nearest neighbor indices and return a pointcloud
neigh::PointCloud<double> neigh::treeKNN::byNumber(int particle,
                                                   size_t nearest) {
  size_t realNeighbors = nearest + 1;
  neigh::PointCloud<double> resultCloud;
  std::vector<size_t> ret_index(realNeighbors);
  std::vector<double> out_dist_sqr(realNeighbors);
  double X = cloud.pts[particle].x;
  double Y = cloud.pts[particle].y;
  double Z = cloud.pts[particle].z;
  const double query_pt[3] = {X, Y, Z};
  // construct a kd-tree index:
  typedef nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor_MD<double, neigh::PointCloud<double>>,
      neigh::PointCloud<double>, 3 /* dim */
      >
      my_tree;
  my_tree index(3 /*dim*/, cloud,
                nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
  index.buildIndex();

#if 0
	// Test resize of dataset and rebuild of index:
	cloud.pts.resize(cloud.pts.size()*0.5);
	index.buildIndex();
#endif

  realNeighbors = index.knnSearch(&query_pt[0], realNeighbors, &ret_index[0],
                                  &out_dist_sqr[0]);

  // In case of less points in the tree than requested:
  ret_index.resize(realNeighbors);
  out_dist_sqr.resize(realNeighbors);

  std::cout << "knnSearch(): num_results=" << realNeighbors << "\n";
  for (size_t i = 0; i < realNeighbors; i++)
    std::cout << "idx[" << i << "]=" << ret_index[i] << " dist[" << i
              << "]=" << out_dist_sqr[i] << std::endl;
  std::cout << "\n";

  std::cout << "point coordinates"
            << "\n";
  for (size_t i = 0; i < 3; i++)
    std::cout << " coord[" << i << "]=" << cloud.kdtree_get_pt(ret_index[0], i);
  std::cout << "\n";

  // Prepare output
  resultCloud.pts.resize(nearest);
  double resX, resY, resZ;
  for (int i = 0; i < nearest; i++) {
    resX = cloud.kdtree_get_pt(ret_index[i + 1], 0);
    resY = cloud.kdtree_get_pt(ret_index[i + 1], 1);
    resZ = cloud.kdtree_get_pt(ret_index[i + 1], 2);
    resultCloud.pts[i].x = resX;
    resultCloud.pts[i].y = resY;
    resultCloud.pts[i].z = resZ;
  }

  // Prep box for output
  resultCloud.box[0] = frame->parameter->boxx;
  resultCloud.box[1] = frame->parameter->boxy;
  resultCloud.box[2] = frame->parameter->boxz;
  // Send resultCloud out
  return resultCloud;
}

// Helper
bool neigh::treeKNN::isThere(int iatom, CMolecularSystem *frame) {
  // TODO: Migrate to CGeneric
  double coordX = frame->molecules[iatom].get_posx();
  double coordY = frame->molecules[iatom].get_posy();
  double coordZ = frame->molecules[iatom].get_posz();
  std::array<double, 3> coord = {coordX, coordY, coordZ};
  for (int i = 0; i < 3; i++) {
    if (coordHigh[i] == coordLow[i]) {
      return true;
    } else if (coord[i] >= coordLow[i] && coord[i] <= coordHigh[i]) {
      return true;
    } else {
      return false;
    }
  }
  return false;
}

// Dummy driver (int for error handling)
// Since the main function reads the traj file, nop is known
int neigh::treeKNN::initKNN(int nop, std::string filename, int frameNum,
                            int typeI) {
  // nsteps must be greater than the frame number here (single frame handling)
  int nstep = frameNum + 1;
  // Prepares the frame
  this->prepFrame(nop, filename);
  // Use nsteps (dummy)
  this->frame->parameter->nsteps = nstep;
  // Read the file
  this->frame->readParticleFile(frameNum);
  // Now we populate the cloud (make int and do error handling)
  this->populateCloud(typeI);
  return 1;
}
