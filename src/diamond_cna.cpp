//-----------------------------------------------------------------------------------
// d-SEAMS
// Copyright (c) d-SEAMS developers
//-----------------------------------------------------------------------------------

#include <diamond_cna.hpp>
#include <local_symmetry.hpp>

#include <vector>

namespace symm {

void diamondFromCsr(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                    const std::vector<int> &offsets, const std::vector<int> &cols,
                    std::vector<int> &kind) {
  const int n = yCloud.nop;
  kind.assign(static_cast<std::size_t>(std::max(n, 0)), -1);
  if (n <= 0 || yCloud.box.size() < 3) {
    return;
  }
  // The common-neighbour walk on this CSR was the slow diamond test.
  // The shell kernel is the per-atom label, including an ideal shell
  // that does not build the signature. The CSR is the caller's cutoff
  // list; the kernel builds the same 3.5 A bonds.
  (void)offsets;
  (void)cols;
  std::vector<double> x(static_cast<std::size_t>(n));
  std::vector<double> y(static_cast<std::size_t>(n));
  std::vector<double> z(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(i)];
    x[static_cast<std::size_t>(i)] = p.x;
    y[static_cast<std::size_t>(i)] = p.y;
    z[static_cast<std::size_t>(i)] = p.z;
  }
  double lat[3][3] = {};
  lat[0][0] = yCloud.box[0];
  lat[1][1] = yCloud.box[1];
  lat[2][2] = yCloud.box[2];
  if (yCloud.box.size() >= 6) {
    lat[1][0] = yCloud.box[3];
    lat[2][0] = yCloud.box[4];
    lat[2][1] = yCloud.box[5];
  }
  classifyIceArrays(x.data(), y.data(), z.data(), n, lat, kind.data());
}

} // namespace symm
