// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
// One CHILL+ census: hexagonal count, primitive six-rings, mean local q6.

#pragma once

#include "seams.h"

#include <vector>

struct FrameCensus {
  int hexagonal = 0;
  int six = 0;
  int cubic = 0;
  double q6 = 0.0;
  bool ok = false;
};

inline FrameCensus frame_census(const double *xyz, int n, const double *box9, double cutoff) {
  FrameCensus out;
  if (xyz == nullptr || box9 == nullptr || n <= 0) {
    return out;
  }
  seams_system *sys = seams_system_create();
  if (sys == nullptr ||
      seams_system_set_frame(sys, xyz, n, box9, nullptr, nullptr) != SEAMS_OK) {
    seams_system_destroy(sys);
    return out;
  }
  std::vector<int> labels(static_cast<std::size_t>(n));
  int counts[9] = {};
  std::vector<double> q6(static_cast<std::size_t>(n));
  if (seams_chill_plus(sys, cutoff, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels.data(),
                       counts) != SEAMS_OK ||
      seams_q4q6(sys, cutoff, 0, nullptr, q6.data()) != SEAMS_OK) {
    seams_system_destroy(sys);
    return out;
  }
  double sum = 0.0;
  for (double value : q6) {
    sum += value;
  }
  seams_rings *rings = nullptr;
  if (seams_rings_find(sys, cutoff, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) != SEAMS_OK) {
    seams_system_destroy(sys);
    return out;
  }
  int hist[7] = {};
  if (seams_rings_histogram(rings, hist, 7) != SEAMS_OK) {
    seams_rings_destroy(rings);
    seams_system_destroy(sys);
    return out;
  }
  seams_rings_destroy(rings);
  seams_system_destroy(sys);
  out.hexagonal = counts[SEAMS_LABEL_HEXAGONAL];
  out.cubic = counts[SEAMS_LABEL_CUBIC];
  out.six = hist[6];
  out.q6 = sum / static_cast<double>(n);
  out.ok = true;
  return out;
}
