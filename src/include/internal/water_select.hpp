// SPDX-License-Identifier: MIT
#ifndef SEAMS_WATER_SELECT_H_
#define SEAMS_WATER_SELECT_H_

#include <mol_sys.hpp>

#include <string>
#include <vector>

/** Water-oxygen selection for a mixed frame.
 *
 * Neighbours, CHILL, CHILL+, rings and cages run on the oxygen atoms.
 * Ions and hydrogens stay out of that graph. A requested type or molecule
 * that is not in the frame is an error, as is a selection that keeps nothing.
 * The messages are the strings the tests match.
 */

namespace water {

enum class MolMode { None, Ids, Water };

struct Request {
  int type = 0;
  std::vector<int> types;
  std::vector<int> molecules;
  MolMode mol_mode = MolMode::None;
  std::vector<int> ion_types;
};

struct Selection {
  std::vector<char> keep;
  std::vector<int> index;
  std::vector<int> oxygen_types;
  int n_oxygen = 0;
};

struct IonSite {
  double x = 0;
  double y = 0;
  double z = 0;
  int type = 0;
};

struct IonReport {
  int ice = 0;
  int liquid = 0;
  double front = -1;
  bool any = false;
};

bool hasMultiAtomMolecule(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud);

Selection selectOxygens(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud,
    const Request &request);

molSys::PointCloud<molSys::Point<double>, double> subsetCloud(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud,
    const Selection &selection);

/** Historical dump guess: type 2 when any atom has it, otherwise type 1.
 *  A miss leaves an empty cloud and keeps the box. */
molSys::PointCloud<molSys::Point<double>, double> keepHistoricalOxygen(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud);

std::vector<IonSite> collectIons(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud,
    const Selection &selection);

/** oxygen is the labelled water subset. front is the nearest ice oxygen,
 *  in the box units, or -1 when no ion sits next to ice. */
IonReport ionFront(
    const molSys::PointCloud<molSys::Point<double>, double> &oxygen,
    const std::vector<IonSite> &ions);

bool iceLabel(molSys::atom_state_type label);

} // namespace water

#endif
