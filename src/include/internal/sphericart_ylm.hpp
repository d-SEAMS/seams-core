#ifndef SEAMS_SPHERICART_YLM_H_
#define SEAMS_SPHERICART_YLM_H_

// Sphericart (Bigi, Fraux, Browning, Ceriotti, J. Chem. Phys. 159,
// 064802 (2023); 10.1063/5.0156307) evaluates real Y_lm from
// Cartesian coordinates. Host Steinhardt pass 1 can batch every bond.
// The OpenMP-target kernel keeps the closed-form ylmAll.

namespace seams {
namespace sphericart_ylm {

bool available();

// xyz is nVec * 3 Cartesian components (need not be unit). ylmOut is
// nVec * (2*orderL+1) * 2 interleaved complex values in the same m
// order as ylmAll (m = -l .. +l). Real spherical harmonics are packed
// into the real parts; imaginary parts are zero. Returns 0 on success.
int ylmCartesian(int orderL, const double *xyz, int nVec, double *ylmOut);

} // namespace sphericart_ylm
} // namespace seams

#endif
