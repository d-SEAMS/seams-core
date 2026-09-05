#include <seams_c_api.h>

#include <bop.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <vector>

int seams_chill_plus(const double *xyz, int n, const double *box, int *labels) {
  if (xyz == nullptr || box == nullptr || labels == nullptr || n <= 0) {
    return 1;
  }
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.nop = n;
  cloud.currentFrame = 1;
  cloud.box = {box[3], box[4], box[5]};
  cloud.boxLow = {box[0], box[1], box[2]};
  cloud.pts.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    auto &p = cloud.pts[static_cast<std::size_t>(i)];
    p.x = xyz[3 * i];
    p.y = xyz[3 * i + 1];
    p.z = xyz[3 * i + 2];
    p.type = 1;
    p.atomID = i + 1;
    p.molID = i + 1;
    cloud.idIndexMap[i + 1] = i;
  }
  const auto nList = nneigh::kNearestNeighbourList(cloud, 4, 5.0, 1, true);
  chill::getCorrelPlus(cloud, nList, false);
  chill::getIceTypePlusNoPrint(cloud, nList, false);
  for (int i = 0; i < n; i++) {
    labels[i] = static_cast<int>(cloud.pts[static_cast<std::size_t>(i)].iceType);
  }
  return 0;
}
