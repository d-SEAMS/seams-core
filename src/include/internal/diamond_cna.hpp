#ifndef SEAMS_DIAMOND_CNA_H_
#define SEAMS_DIAMOND_CNA_H_

// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT

#include <mol_sys.hpp>

#include <vector>

// Per-atom ice label from the shell kernel. A clean cubic shell is 1
// (ice Ic), a clean hexagonal shell is 0 (ice Ih), stacking-disordered
// ice is 2, and anything else is -1. The cutoff CSR is accepted so older
// callers keep their signature; the kernel builds the 3.5 A shell itself.

namespace symm {

void diamondFromCsr(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                    const std::vector<int> &offsets, const std::vector<int> &cols,
                    std::vector<int> &kind);

} // namespace symm

#endif
