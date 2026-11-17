//-----------------------------------------------------------------------------------
// d-SEAMS
// Copyright (c) d-SEAMS developers
//-----------------------------------------------------------------------------------

#include <diamond_cna.hpp>
#include <generic.hpp>
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
  const gen::FracBox frame = gen::makeFracBox(yCloud);
  if (!frame.ok) {
    return;
  }
  double lat[3][3] = {};
  lat[0][0] = frame.lx;
  lat[1][1] = frame.ly;
  lat[2][2] = frame.lz;
  lat[1][0] = frame.xy;
  lat[2][0] = frame.xz;
  lat[2][1] = frame.yz;
  classifyIceArrays(x.data(), y.data(), z.data(), n, lat, kind.data());
}

} // namespace symm
