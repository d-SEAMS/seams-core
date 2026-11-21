/* Copyright (c) d-SEAMS developers
 * seams::unwrap_positions and seams::unwrap_centres must refuse an xyz span
 * that is not nframes * 3 * natoms long, and a count whose product
 * overflows size_t, before any allocation or C call.
 */
#include "seams.hpp"

#include <cstddef>
#include <cstdio>
#include <exception>
#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace {

int failures = 0;

void expect_range(const char *label, const std::function<void()> &call, const char *needle) {
  try {
    call();
  } catch (const seams::Error &err) {
    if (err.code() != SEAMS_ERR_RANGE || std::string(err.what()).find(needle) == std::string::npos) {
      std::fprintf(stderr, "FAIL %s: status %d, message '%s'\n", label, err.code(), err.what());
      ++failures;
    } else {
      std::printf("ok %s: %s\n", label, err.what());
    }
    return;
  } catch (const std::exception &err) {
    std::fprintf(stderr, "FAIL %s: threw %s, not seams::Error\n", label, err.what());
    ++failures;
    return;
  }
  std::fprintf(stderr, "FAIL %s: accepted\n", label);
  ++failures;
}

}  // namespace

int main() {
  const std::vector<double> box9 = {0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 0.0, 0.0, 0.0};
  const std::size_t nframes = 3;
  const std::size_t natoms = 2;
  std::vector<double> xyz(nframes * 3 * natoms);
  for (std::size_t i = 0; i < xyz.size(); ++i) {
    xyz[i] = 1.0 + 0.1 * static_cast<double>(i);
  }
  const std::vector<int> mols = {1, 1};

  const auto good = seams::unwrap_positions(xyz, nframes, natoms, box9);
  if (good.size() != xyz.size()) {
    std::fprintf(stderr, "FAIL exact span: %zu values\n", good.size());
    ++failures;
  }
  const auto good_c = seams::unwrap_centres(xyz, nframes, natoms, box9, mols);
  if (good_c.mol_ids.size() != 1 || good_c.centres.size() != nframes * 3) {
    std::fprintf(stderr, "FAIL exact span centres: %zu ids, %zu values\n", good_c.mol_ids.size(),
                 good_c.centres.size());
    ++failures;
  }

  // Short spans are a heap read past the end unless refused.
  const std::vector<double> short_xyz(xyz.begin(), xyz.end() - 1);
  const std::vector<double> long_xyz = [&] {
    auto v = xyz;
    v.push_back(0.5);
    return v;
  }();
  const std::size_t most = std::numeric_limits<std::size_t>::max();
  // (most / 3 + 2) * 3 wraps to 5; most / 2 * 3 wraps to a length no vector holds.
  const std::size_t wraps_small = most / 3 + 2;
  const std::size_t wraps_huge = most / 2;
  const char *length = "xyz length is not nframes * 3 * natoms";
  const char *overflow = "nframes * 3 * natoms overflows size_t";

  expect_range("unwrap_positions short span",
               [&] { (void)seams::unwrap_positions(short_xyz, nframes, natoms, box9); }, length);
  expect_range("unwrap_positions long span",
               [&] { (void)seams::unwrap_positions(long_xyz, nframes, natoms, box9); }, length);
  expect_range("unwrap_positions empty span",
               [&] { (void)seams::unwrap_positions({}, nframes, natoms, box9); }, length);
  expect_range("unwrap_positions nframes wraps small",
               [&] { (void)seams::unwrap_positions(xyz, wraps_small, 1, box9); }, overflow);
  expect_range("unwrap_positions nframes wraps huge",
               [&] { (void)seams::unwrap_positions(xyz, wraps_huge, 1, box9); }, overflow);
  expect_range("unwrap_positions natoms over SIZE_MAX / 3",
               [&] { (void)seams::unwrap_positions(xyz, 1, most / 3 + 1, box9); }, overflow);

  expect_range("unwrap_centres short span",
               [&] { (void)seams::unwrap_centres(short_xyz, nframes, natoms, box9, mols); }, length);
  expect_range("unwrap_centres long span",
               [&] { (void)seams::unwrap_centres(long_xyz, nframes, natoms, box9, mols); }, length);
  expect_range("unwrap_centres nframes wraps small",
               [&] { (void)seams::unwrap_centres(xyz, wraps_small, 2, box9, mols); }, overflow);
  expect_range("unwrap_centres nframes wraps huge",
               [&] { (void)seams::unwrap_centres(xyz, wraps_huge, 2, box9, mols); }, overflow);

  if (failures != 0) {
    std::fprintf(stderr, "%d failures\n", failures);
    return 1;
  }
  std::printf("seams.hpp trajectory bounds ok\n");
  return 0;
}
