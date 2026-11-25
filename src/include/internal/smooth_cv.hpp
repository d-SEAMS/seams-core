//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#ifndef SEAMS_SMOOTH_CV_H_
#define SEAMS_SMOOTH_CV_H_

#include <array>
#include <cstddef>
#include <vector>

#include <mol_sys.hpp>

/** @file smooth_cv.hpp
 *  @brief Differentiable ring, cage and gyration collective variables.
 */

/**
 *  @addtogroup smooth
 *  @{
 */

/** @brief Smooth collective variables on the primitive-ring network.
 *  @details Every bond of the graph carries a weight s(r) from a rational
 *  switching function with a C2 taper to zero at d_max. Every simple cycle
 *  of at most maxRing members found by primitive::countAllRingsFromIndex on
 *  the graph at d_max + skin is a candidate ring. Its weight is the product
 *  of its bond weights times one factor (1 - W_P) for every shortcut path P,
 *  a simple path between two ring members that is shorter than their
 *  separation around the ring, where W_P is the product of the path's bond
 *  weights. With weights 0 and 1 the rings of weight one are exactly the
 *  shortest-path rings of primitive::ringNetwork.
 *
 *  Every ring weight, ring participation and hexagon-pair term is a fixed
 *  polynomial in the pair weights of the complete graph. A pair outside the
 *  neighbour list has weight zero, and a term that contains a zero bond
 *  weight, or a factor (1 - 0), contributes what it would contribute if the
 *  pair were listed. Those values therefore do not depend on which pairs
 *  beyond d_max the list holds, and they are as smooth as s(r) at every
 *  geometry, including where a bond crosses the list cutoff or a ring
 *  changes primitivity. The gyration values also depend on the unwrapped
 *  positions: they are as smooth only for a compact selection of positive
 *  total weight, and the volume is not differentiable where det G = 0.
 */
namespace smooth {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

/** @brief PLUMED RATIONAL switching function times a C2 taper.
 *  @details For x = (r - d0) / r0 the rational part is
 *  (1 - x^n) / (1 - x^m), evaluated as the ratio of the geometric sums
 *  1 + x + ... + x^(n-1) and 1 + x + ... + x^(m-1), which has no removable
 *  singularity at x = 1, and above x = 1 as y^(m-n) times the same ratio
 *  in y = 1 / x. It is 1 for r at or below d0. The taper is 1 below
 *  dOn, 0 at and beyond dMax, and 1 - u^3 (10 - 15 u + 6 u^2) with
 *  u = (r - dOn) / (dMax - dOn) in between, so the product and its first
 *  two derivatives vanish at dMax. Near d0 the rational part is
 *  1 - x^n + O(x^m), so s is C^(n-1) and not C^n at d0 > 0; at d0 = 0
 *  and n = 1, s(|v|) has a cone at v = 0.
 */
struct RationalSwitch {
  double r0 = 1.0;
  double d0 = 0.0;
  int n = 6;
  int m = 12;
  double dOn = 0.0;
  double dMax = 1.0;

  //! Throws std::invalid_argument unless r0 > 0, n >= 1, m > n,
  //! 0 <= d0 <= dOn < dMax and every value is finite. m == 0 is read as 2n
  //! by makeSwitch before this check.
  void validate() const;

  //! s(r), finite and not negative; dsdr receives ds/dr, the right derivative
  //! at r = d0. No power of x above 1 is formed, so only dsdr can leave
  //! double range, when r0 or dMax - dOn is near the smallest doubles.
  double operator()(double r, double &dsdr) const;
};

//! Builds and validates a switch. m == 0 selects m = 2n, as in PLUMED,
//! and throws std::invalid_argument when 2n does not fit in an int.
RationalSwitch makeSwitch(double r0, double d0, int n, int m, double dOn,
                          double dMax);

/** @brief Candidate rings, shortcut paths and hexagon pairs of one graph.
 *  @details Fixed between neighbour-list rebuilds. Edge e joins
 *  edges[e][0] < edges[e][1]. Ring k lists its members in
 *  ringAtoms[ringOff[k] .. ringOff[k + 1]], in cycle order from its lowest
 *  index toward the lesser of that member's ring neighbours, and the edge
 *  from member t to member t + 1 (cyclically) in ringEdges at the same
 *  offsets. The shortcut paths of ring k are paths pathOff[k] ..
 *  pathOff[k + 1]; path p is the edge list
 *  pathEdges[pathEdgeOff[p] .. pathEdgeOff[p + 1]].
 */
struct Topology {
  int nAtoms = 0;
  int maxRing = 0;
  std::vector<std::array<int, 2>> edges;
  std::vector<int> adjOff;
  std::vector<int> adjAtom;
  std::vector<int> adjEdge;
  std::vector<int> ringOff{0};
  std::vector<int> ringAtoms;
  std::vector<int> ringEdges;
  std::vector<int> pathOff{0};
  std::vector<int> pathEdgeOff{0};
  std::vector<int> pathEdges;

  /** @brief Two vertex-disjoint candidate hexagons a < b that may pass the
   *  HC basal test. bond[6 * i + j] is the edge between member i of a and
   *  member j of b, or -1 when the pair is not in the graph. */
  struct HexPair {
    int a = 0;
    int b = 0;
    std::array<int, 36> bond{};
  };
  std::vector<HexPair> hexPairs;

  [[nodiscard]] int ringCount() const { return static_cast<int>(ringOff.size()) - 1; }
  [[nodiscard]] int ringSize(int k) const {
    return ringOff[static_cast<std::size_t>(k) + 1] - ringOff[static_cast<std::size_t>(k)];
  }
  //! Edge id of the pair (i, j), or -1
  [[nodiscard]] int edgeOf(int i, int j) const;
};

//! Topology of a graph given as rows by index, each row starting with the
//! vertex itself (the layout of primitive::ringNetwork). The rows must be
//! symmetric. Throws std::invalid_argument for maxRing outside 3 .. 10.
Topology buildTopology(const std::vector<std::vector<int>> &rows, int maxRing);

//! Ring weights for bond weights s, one per edge
std::vector<double> ringWeights(const Topology &topo, const std::vector<double> &s);

//! Adds d(sum_k wbar[k] w[k]) / ds to sbar
void ringWeightsAdjoint(const Topology &topo, const std::vector<double> &s,
                        const std::vector<double> &wbar, std::vector<double> &sbar);

//! Soft count of hexagon pairs that pass the HC basal test, from bond
//! weights s and ring weights w. With sbar and wbar non-null the adjoint of
//! the count is added to them.
double hexBasalPairs(const Topology &topo, const std::vector<double> &s,
                     const std::vector<double> &w, std::vector<double> *sbar,
                     std::vector<double> *wbar);

//! Soft ring participation of every atom: 1 - prod (1 - w_k) over the
//! candidate rings through the atom whose size is size (0 is every size)
std::vector<double> ringParticipation(const Topology &topo, const std::vector<double> &w,
                                      int size);

//! Adds d(sum_i pbar[i] p[i]) / dw to wbar
void ringParticipationAdjoint(const Topology &topo, const std::vector<double> &w, int size,
                              const std::vector<double> &pbar, std::vector<double> &wbar);

/** @brief The collective variables of one frame, with their gradients.
 *  @details The first three are as smooth as the switch; see the namespace
 *  notes for the gyration kinds. A state is as smooth as its clauses. */
enum class Kind {
  RingCount = 0,
  RingAtoms = 1,
  HexBasalPairs = 2,
  GyrationRg2 = 3,
  GyrationVolume = 4,
  State = 5
};

/** @brief Switch, candidate topology, rebuild reference and state clauses.
 *  @details rebuild() lists the selected atoms within dMax + skin with
 *  linkcell and fixes the candidate topology. value() equals the value
 *  after a fresh rebuild (for the gyration kinds, of a compact selection)
 *  while no selected atom has moved more than skin / 2 from its position
 *  at the rebuild and the cell is unchanged; stale() reports when that
 *  fails.
 */
class Workspace {
public:
  Workspace(const RationalSwitch &sw, int maxRing);

  //! Lists the graph and builds the topology. selected holds frame
  //! indices, ascending and unique. Throws std::invalid_argument when the
  //! cell has no volume or dMax + skin is not below half its narrowest width.
  void rebuild(const Cloud &cloud, std::vector<int> selected, double skin);

  [[nodiscard]] bool built() const { return built_; }
  [[nodiscard]] bool stale(const Cloud &cloud) const;
  [[nodiscard]] const Topology &topology() const { return topo_; }
  [[nodiscard]] const RationalSwitch &switching() const { return sw_; }
  [[nodiscard]] int maxRing() const { return maxRing_; }
  [[nodiscard]] int stateCount() const { return static_cast<int>(states_.size()); }

  //! Throws std::invalid_argument when kind and arg do not name a variable
  //! of this workspace.
  void checkKind(Kind kind, int arg) const;

  //! Value of one variable. With grad non-null, 3 * cloud.nop entries of
  //! d value / d x are written, zero for unselected atoms. Throws
  //! std::range_error, writing nothing, when the value or a gradient entry
  //! is not finite.
  double value(const Cloud &cloud, Kind kind, int arg, double *grad);

  //! Appends a factor to state `state` (state == stateCount() opens a new
  //! state). The factor is sigma((v - threshold) / width) when above is
  //! true and 1 - sigma(...) otherwise, sigma the logistic function. With
  //! kind State and width 0 the factor is chi_arg or 1 - chi_arg.
  void addClause(int state, Kind kind, int arg, double threshold, double width, bool above);

  //! chi of every state; with grad non-null, stateCount() rows of
  //! 3 * cloud.nop entries. Throws std::range_error, writing nothing, when
  //! a chi or a gradient entry is not finite.
  void states(const Cloud &cloud, double *chi, double *grad);

private:
  struct Clause {
    Kind kind;
    int arg;
    double threshold;
    double width;
    bool above;
  };

  void edgeWeights(const Cloud &cloud);
  double channel(const Cloud &cloud, Kind kind, int arg, std::vector<double> *grad);
  double gyration(const Cloud &cloud, bool volume, int arg, std::vector<double> *grad);
  void stateValues(const Cloud &cloud, std::vector<double> &chi,
                   std::vector<std::vector<double>> *grads);
  void edgesToAtoms(const std::vector<double> &sbar, std::vector<double> &grad) const;

  RationalSwitch sw_;
  int maxRing_ = 6;
  bool built_ = false;
  Topology topo_;
  std::vector<int> selected_;
  std::vector<double> reference_;
  std::vector<double> box_;
  std::vector<double> boxLow_;
  double skin_ = 0.0;
  std::vector<std::vector<Clause>> states_;
  std::vector<double> s_;
  std::vector<double> ds_;
  std::vector<std::array<double, 3>> unit_;
};

} // namespace smooth

/** @} */
#endif // SEAMS_SMOOTH_CV_H_
