#ifndef SEAMS_CHILL_OFFLOAD_H_
#define SEAMS_CHILL_OFFLOAD_H_

#include <mol_sys.hpp>
#include <vector>

/** @file chill_offload.hpp
 *  @brief CHILL+ labels: host classifyBonds or flattened OpenMP target.
 *
 *  SEAMS_OFFLOAD=0 uses getCorrelPlus + getIceTypePlusNoPrint.
 *  Any other value uses the flattened four-nearest CHILL+ tables, on
 *  the OpenMP target when a device exists.
 */

namespace chill {

struct ChillPlusResult {
  std::vector<int> iceType;
  bool usedDevice = false;
};

[[nodiscard]] bool preferOffload();

/** Host CHILL+ (classifyBonds + getIceTypePlusNoPrint). */
[[nodiscard]] ChillPlusResult hostChillPlus(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList);

/** Flattened CHILL+. Target device when compiled with offload. */
[[nodiscard]] ChillPlusResult specializedChillPlus(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList);

[[nodiscard]] ChillPlusResult chillPlus(
    molSys::PointCloud<molSys::Point<double>, double> &yCloud,
    const std::vector<std::vector<int>> &nList);

} // namespace chill

#endif
