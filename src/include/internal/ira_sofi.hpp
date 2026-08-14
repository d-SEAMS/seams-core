#ifndef SEAMS_IRA_SOFI_H_
#define SEAMS_IRA_SOFI_H_

#include <Eigen/Core>
#include <string>
#include <vector>

// IRA (Gunde, Salles, Hemeryck, Martin-Samos, JCIM 2021,
// 10.1021/acs.jcim.1c00567) solves rotation + translation + permutation.
// SOFI (Gunde, Salles, Grisanti, Hemeryck, Martin-Samos, JCP 2024,
// 10.1063/5.0215689) returns the degenerate self-matches as the point group.
// Both live in libira. A build without the library still compiles these
// declarations; match and pointGroup then return 1.

namespace ira {

struct Match {
  Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();
  Eigen::Vector3d translation = Eigen::Vector3d::Zero();
  std::vector<int> assignment;
  std::vector<double> quat;
  double hausdorff = 0.0;
  double rmsd = 0.0;
};

bool available();

// Overlay target onto ref. Both are n x 3 (or m x 3). Equal size is the
// cage/prism case; unequal size is a fragment match. Returns 0 on success.
[[nodiscard]] int match(const Eigen::MatrixXd &ref,
                        const Eigen::MatrixXd &target, Match &out,
                        double kmaxFactor = 1.8);

// Point-group operations of a centered n x 3 cloud. Returns 0 on success.
struct PointGroup {
  std::string symbol;
  int nOperations = 0;
  std::vector<Eigen::Matrix3d> operations;
};

[[nodiscard]] int pointGroup(const Eigen::MatrixXd &coords, PointGroup &out,
                             double symThr = 0.05);

// Horn-shaped result for callers that already store quat + rmsd.
// Returns true when IRA produced a match.
bool orient(const Eigen::MatrixXd &ref, const Eigen::MatrixXd &target,
            std::vector<double> &quat, double &rmsd);

} // namespace ira

#endif // SEAMS_IRA_SOFI_H_
