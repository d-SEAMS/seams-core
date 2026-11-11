#ifndef SEAMS_DIAMOND_CNA_H_
#define SEAMS_DIAMOND_CNA_H_

// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT

#include <mol_sys.hpp>

#include <vector>

// Diamond common-neighbour test. A clean cubic shell is 1 and a clean
// hexagonal shell is 0. Anything the test cannot decide is -1. The twelve
// second-shell vectors are read from the cutoff CSR, not from a second search.

namespace symm {

void diamondFromCsr(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                    const std::vector<int> &offsets, const std::vector<int> &cols,
                    std::vector<int> &kind);

} // namespace symm

#endif
