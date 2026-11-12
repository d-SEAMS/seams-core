#include <catch2/catch_test_macros.hpp>

#include <bop.hpp>
#include <chill_plus_scheme.hpp>
#include <franzblau.hpp>
#include <neighbours.hpp>
#include <seams_c_api.h>
#include <seams_input.hpp>
#include <water_select.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace {

void pushPoint(std::vector<double> &xyz, double x, double y, double z) {
  xyz.push_back(x);
  xyz.push_back(y);
  xyz.push_back(z);
}

std::vector<double> tetrahedron() {
  std::vector<double> xyz;
  pushPoint(xyz, 10.0, 10.0, 10.0);
  pushPoint(xyz, 12.8, 10.0, 10.0);
  pushPoint(xyz, 9.0667, 12.424, 10.0);
  pushPoint(xyz, 9.0667, 8.788, 11.616);
  pushPoint(xyz, 9.0667, 8.788, 8.384);
  return xyz;
}

const double kBox9[9] = {0, 0, 0, 30, 30, 30, 0, 0, 0};

int ringCensus(molSys::PointCloud<molSys::Point<double>, double> &cloud,
               double cutoff, int type, std::vector<int> &counts) {
  auto ids = nneigh::neighListO(cutoff, cloud, type);
  auto idx = nneigh::neighbourListByIndex(cloud, ids);
  const auto rings = primitive::ringNetwork(idx, 6);
  counts.assign(7, 0);
  for (const auto &ring : rings) {
    if (ring.size() < counts.size()) {
      counts[ring.size()] += 1;
    }
  }
  return 0;
}

molSys::PointCloud<molSys::Point<double>, double>
cloudOf(const std::vector<double> &xyz, const std::vector<int> &types,
        const std::vector<int> &mols) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const int n = static_cast<int>(types.size());
  cloud.nop = n;
  cloud.currentFrame = 1;
  cloud.box = {30, 30, 30};
  cloud.boxLow = {0, 0, 0};
  cloud.pts.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    auto &point = cloud.pts[static_cast<std::size_t>(i)];
    point.x = xyz[static_cast<std::size_t>(3 * i)];
    point.y = xyz[static_cast<std::size_t>(3 * i + 1)];
    point.z = xyz[static_cast<std::size_t>(3 * i + 2)];
    point.type = types[static_cast<std::size_t>(i)];
    point.molID = mols[static_cast<std::size_t>(i)];
    point.atomID = i + 1;
    cloud.idIndexMap[point.atomID] = i;
  }
  return cloud;
}

} // namespace

TEST_CASE("untyped CHILL+ labels an ion, selection leaves it unclassified",
          "[brine][chill]") {
  auto oxygen = tetrahedron();
  std::vector<double> mixed = oxygen;
  pushPoint(mixed, 12.0, 10.0, 10.0);
  const int nOx = 5;
  const int nAll = 6;
  std::vector<int> oxLabels(static_cast<std::size_t>(nOx), -1);
  std::vector<int> oldLabels(static_cast<std::size_t>(nAll), -1);
  std::vector<int> newLabels(static_cast<std::size_t>(nAll), -1);
  REQUIRE(seams_chill_plus_graph_lammps(oxygen.data(), nOx, kBox9, 3.5, "cutoff",
                                        oxLabels.data()) == 0);
  // The old in-situ path treats every coordinate as water. The ion, 2 A
  // from the centre, is inside the 3.5 A shell and receives a water label.
  REQUIRE(seams_chill_plus_graph_lammps(mixed.data(), nAll, kBox9, 3.5, "cutoff",
                                        oldLabels.data()) == 0);
  REQUIRE(oldLabels[5] != 6);

  std::vector<int> types = {1, 1, 1, 1, 1, 3};
  int ice = -1;
  int liquid = -1;
  double front = 0;
  REQUIRE(seams_chill_plus_species(mixed.data(), nAll, kBox9, types.data(),
                                   nullptr, 1, 0, nullptr, 0, 3.5, "cutoff",
                                   newLabels.data(), &ice, &liquid,
                                   &front) == 0);
  for (int i = 0; i < nOx; ++i) {
    REQUIRE(newLabels[static_cast<std::size_t>(i)] ==
            oxLabels[static_cast<std::size_t>(i)]);
  }
  REQUIRE(newLabels[5] == 6);
  REQUIRE(ice + liquid == 1);
  REQUIRE((ice == 0) == (front < 0.0));

  std::vector<int> degree(static_cast<std::size_t>(nAll), 0);
  REQUIRE(seams_oxygen_cutoff_degree(mixed.data(), nAll, kBox9, types.data(),
                                     nullptr, 1, 0, nullptr, 0, 3.5,
                                     degree.data()) == 0);
  REQUIRE(degree[5] == -1);
  std::vector<int> oxDegree(static_cast<std::size_t>(nOx), -2);
  REQUIRE(seams_oxygen_cutoff_degree(oxygen.data(), nOx, kBox9, nullptr, nullptr,
                                     0, 0, nullptr, 0, 3.5,
                                     oxDegree.data()) == 0);
  for (int i = 0; i < nOx; ++i) {
    REQUIRE(degree[static_cast<std::size_t>(i)] ==
            oxDegree[static_cast<std::size_t>(i)]);
  }
}

TEST_CASE("a mixed frame without a type or molecule is an error", "[brine]") {
  auto mixed = tetrahedron();
  pushPoint(mixed, 12.0, 10.0, 10.0);
  std::vector<int> types = {1, 1, 1, 1, 1, 3};
  int label = 0;
  REQUIRE(seams_chill_plus_species(mixed.data(), 6, kBox9, types.data(), nullptr,
                                   0, 0, nullptr, 0, 3.5, "cutoff", &label,
                                   nullptr, nullptr, nullptr) == 1);
  REQUIRE(std::string(seams_c_error()) ==
          "mixed species need a water type or molecule");
}

TEST_CASE("bad type and molecule selections name the miss", "[brine]") {
  auto xyz = tetrahedron();
  std::vector<int> types(5, 1);
  std::vector<int> mols = {1, 1, 1, 1, 2};
  types.push_back(3);
  mols.push_back(3);
  pushPoint(xyz, 20.0, 10.0, 10.0);
  auto cloud = cloudOf(xyz, types, mols);

  auto requireMessage = [&](const water::Request &request, const char *text) {
    bool threw = false;
    try {
      water::selectOxygens(cloud, request);
    } catch (const std::invalid_argument &ex) {
      threw = true;
      REQUIRE(std::string(ex.what()) == text);
    }
    REQUIRE(threw);
  };
  water::Request missingType;
  missingType.type = 9;
  requireMessage(missingType, "atom type 9 is not in the frame");

  water::Request missingMol;
  missingMol.type = 1;
  missingMol.mol_mode = water::MolMode::Ids;
  missingMol.molecules = {9};
  requireMessage(missingMol, "molecule 9 is not in the frame");

  water::Request ionAsOxygen;
  ionAsOxygen.type = 3;
  ionAsOxygen.ion_types = {3};
  requireMessage(ionAsOxygen, "atom type 3 is an ion, not a water oxygen");

  water::Request waterMol;
  waterMol.mol_mode = water::MolMode::Water;
  // Five monatomic oxygens of type 1 plus one ion: molecules are size 1
  // once the shared mol ids above are... mols 1 is shared by four oxygens.
  // That is a multi-atom molecule whose atoms are all type 1, so inference
  // has no singleton. Ask with an explicit type instead.
  waterMol.type = 1;
  const auto kept = water::selectOxygens(cloud, waterMol);
  REQUIRE(kept.n_oxygen == 5);
  REQUIRE(kept.keep[5] == 0);
}

TEST_CASE("molecule water keeps the TIP4P oxygen and drops Na", "[brine]") {
  std::vector<double> xyz;
  pushPoint(xyz, 10, 10, 10);
  pushPoint(xyz, 10.96, 10, 10);
  pushPoint(xyz, 9.76, 10.93, 10);
  pushPoint(xyz, 14, 10, 10);
  std::vector<int> types = {1, 2, 2, 3};
  std::vector<int> mols = {7, 7, 7, 8};
  auto full = cloudOf(xyz, types, mols);
  water::Request request;
  request.mol_mode = water::MolMode::Water;
  const auto selection = water::selectOxygens(full, request);
  REQUIRE(selection.n_oxygen == 1);
  REQUIRE(selection.oxygen_types.size() == 1);
  REQUIRE(selection.oxygen_types[0] == 1);
  REQUIRE(selection.index[0] == 0);
  const auto ions = water::collectIons(full, selection);
  REQUIRE(ions.size() == 1);
  REQUIRE(ions[0].type == 3);

  int labels[4] = {0, 0, 0, 0};
  const int molKeep = 7;
  REQUIRE(seams_chill_plus_species(xyz.data(), 4, kBox9, types.data(), mols.data(),
                                   0, 2, &molKeep, 1, 3.5, "cutoff", labels,
                                   nullptr, nullptr, nullptr) == 0);
  REQUIRE(labels[1] == 6);
  REQUIRE(labels[2] == 6);
  REQUIRE(labels[3] == 6);
}

TEST_CASE("inferred oxygens match the type-2 read on exampleTraj", "[brine]") {
  molSys::PointCloud<molSys::Point<double>, double> filtered;
  molSys::PointCloud<molSys::Point<double>, double> all;
  filtered = sinp::readLammpsTrjO("traj/exampleTraj.lammpstrj", 1, filtered, 2);
  all = sinp::readLammpsTrj("traj/exampleTraj.lammpstrj", 1, all);
  REQUIRE(filtered.nop > 0);
  REQUIRE(all.nop > filtered.nop);
  water::Request request;
  request.mol_mode = water::MolMode::Water;
  const auto selection = water::selectOxygens(all, request);
  auto oxygen = water::subsetCloud(all, selection);
  REQUIRE(oxygen.nop == filtered.nop);
  for (int i = 0; i < oxygen.nop; ++i) {
    REQUIRE(oxygen.pts[static_cast<std::size_t>(i)].atomID ==
            filtered.pts[static_cast<std::size_t>(i)].atomID);
  }
  auto listA = nneigh::neighListO(3.5, filtered, 2);
  auto listB = nneigh::neighListO(3.5, oxygen, 2);
  REQUIRE(listA == listB);
  chill::applyChillPlus(filtered, listA, chill::chillPlusRule(""), false);
  chill::applyChillPlus(oxygen, listB, chill::chillPlusRule(""), false);
  for (int i = 0; i < oxygen.nop; ++i) {
    REQUIRE(oxygen.pts[static_cast<std::size_t>(i)].iceType ==
            filtered.pts[static_cast<std::size_t>(i)].iceType);
  }
  const auto qA = chill::steinhardtQl(filtered, listA, 6);
  const auto qB = chill::steinhardtQl(oxygen, listB, 6);
  REQUIRE(qA.ql.size() == qB.ql.size());
  double maxAbs = 0.0;
  for (std::size_t i = 0; i < qA.ql.size(); ++i) {
    maxAbs = std::max(maxAbs, std::fabs(qA.ql[i] - qB.ql[i]));
    maxAbs = std::max(maxAbs, std::fabs(qA.qlBar[i] - qB.qlBar[i]));
  }
  REQUIRE(maxAbs <= 1e-12);
  std::vector<int> ringsA;
  std::vector<int> ringsB;
  ringCensus(filtered, 3.5, 2, ringsA);
  ringCensus(oxygen, 3.5, 2, ringsB);
  REQUIRE(ringsA == ringsB);
}

TEST_CASE("a single-species cloud matches the untyped LAMMPS call", "[brine]") {
  auto xyz = tetrahedron();
  const int n = 5;
  std::vector<int> types(static_cast<std::size_t>(n), 1);
  std::vector<int> plain(static_cast<std::size_t>(n), -1);
  std::vector<int> typed(static_cast<std::size_t>(n), -2);
  REQUIRE(seams_chill_plus_graph_lammps(xyz.data(), n, kBox9, 3.5, "cutoff",
                                        plain.data()) == 0);
  REQUIRE(seams_chill_plus_species(xyz.data(), n, kBox9, types.data(), nullptr, 0,
                                   0, nullptr, 0, 3.5, "cutoff", typed.data(),
                                   nullptr, nullptr, nullptr) == 0);
  REQUIRE(plain == typed);
  REQUIRE(seams_chill_plus_species(xyz.data(), n, kBox9, nullptr, nullptr, 0, 0,
                                   nullptr, 0, 3.5, "mutual", typed.data(),
                                   nullptr, nullptr, nullptr) == 0);
  std::vector<int> mutual(static_cast<std::size_t>(n), -3);
  REQUIRE(seams_chill_plus_graph_lammps(xyz.data(), n, kBox9, 3.5, "mutual",
                                        mutual.data()) == 0);
  REQUIRE(typed == mutual);
}

// Two oxygens from the 1 mol/kg, 1000-water brine, each with the two-hop
// shell unwrapped into a large box. Nguyen and Molinero, J. Phys. Chem. B
// 2015, 119, 9369, Table 1: interfacial clathrate is three eclipsed bonds,
// any number of staggered bonds, and four neighbours inside 3.5 A. Table 2:
// eclipsed is -0.35 <= c(i,j) <= 0.25 and staggered is c(i,j) <= -0.8.
// The first oxygen has c = -0.055, 0.234, -0.113, 0.750, so E=3 and one
// bond out of range. That is interfacial clathrate. mdapy 2.0 calls it
// other. The second has one eclipsed bond and one staggered bond. That is
// liquid. mdapy 2.0 calls it interfacial hydrate.
void pushShifted(std::vector<double> &xyz, double x, double y, double z) {
  xyz.push_back(x + 20.0);
  xyz.push_back(y + 20.0);
  xyz.push_back(z + 20.0);
}

int labelAt(const std::vector<double> &xyz, int index) {
  const int n = static_cast<int>(xyz.size() / 3);
  std::vector<int> types(static_cast<std::size_t>(n), 1);
  std::vector<int> labels(static_cast<std::size_t>(n), -1);
  const double box9[9] = {0, 0, 0, 50, 50, 50, 0, 0, 0};
  REQUIRE(seams_chill_plus_species(xyz.data(), n, box9, types.data(), nullptr, 1,
                                   0, nullptr, 0, 3.5, "cutoff", labels.data(),
                                   nullptr, nullptr, nullptr) == 0);
  return labels[static_cast<std::size_t>(index)];
}

TEST_CASE("three eclipsed bonds are interfacial clathrate, one is liquid",
          "[brine][chill][paper]") {
  std::vector<double> threeEclipsed;
  const double clusterE[][3] = {
      {-0.40390000, -3.28911000, -3.00309000},
      {2.63560000, -0.14011000, -3.01262000},
      {-3.42890000, -3.17064000, -0.24502000},
      {-0.17350000, -3.41073000, -0.01629000},
      {0.01420000, -3.59929000, 2.98382000},
      {0.00000000, 0.00000000, 0.00000000},
      {-0.25480000, 0.40561000, 2.74172000},
      {2.62380000, -3.34077000, 0.07012000},
      {3.01770000, 0.09883000, 0.09341000},
      {2.95420000, 0.13878000, 3.04742000},
      {5.92800000, -0.00068000, -0.00017000},
      {-3.84320000, 3.67924000, 0.25871000},
      {-0.48080000, 3.22815000, 0.40350000},
      {-0.21090000, 3.53271000, 3.19752000},
      {-0.81990000, 6.03011000, -0.54107000},
      {2.87800000, 3.22007000, -0.11282000},
      {-0.25460000, -6.22597600, -0.05595000},
  };
  for (const auto &p : clusterE) {
    pushShifted(threeEclipsed, p[0], p[1], p[2]);
  }
  // seams label 5 is interClathrate.
  REQUIRE(labelAt(threeEclipsed, 5) == 5);

  std::vector<double> oneEclipsed;
  const double clusterL[][3] = {
      {0.19365000, -6.01050000, 0.78753000},
      {-0.03092000, -2.98450000, -2.43420000},
      {0.29834000, -2.75870000, 0.84058000},
      {0.08895000, -0.13310000, -2.66461000},
      {0.00000000, 0.00000000, 0.00000000},
      {3.18243000, -3.04970000, 0.49707000},
      {3.23489000, 0.32100000, -2.39453000},
      {-6.00450600, 0.18880000, 0.73629000},
      {-2.99860600, -2.81700000, 0.68921000},
      {-2.92910600, 0.21910000, -2.51959000},
      {-3.29520600, 0.28500000, 0.64232000},
      {0.05733000, 3.44860000, -2.40260000},
      {0.15900000, 3.38230000, 0.49181000},
      {0.06704000, 6.37860000, 0.58690000},
      {3.23097000, 3.44040000, 0.72609000},
      {-3.17410600, 3.32980000, 0.67909000},
      {0.02458000, -3.04620000, 3.84724000},
      {-2.92640600, 0.31480000, 3.70024000},
      {0.17616000, 3.46690000, 3.97215000},
      {0.41108000, 0.05340000, -5.91502600},
  };
  for (const auto &p : clusterL) {
    pushShifted(oneEclipsed, p[0], p[1], p[2]);
  }
  // seams label 2 is water. One eclipsed bond is not Table 1's 3E row.
  REQUIRE(labelAt(oneEclipsed, 4) == 2);
}
