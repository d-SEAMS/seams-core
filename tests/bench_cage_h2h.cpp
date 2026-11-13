/*
** Time the cage search that starts from a primitive ring network:
** hexagonal cages, double-diamond cages, and the clathrate signatures
** 5^12, 5^12 6^2 and 5^12 6^4. Ring construction is outside the timer.
*/
#include <cage.hpp>
#include <cage_affiliation.hpp>
#include <cage_enum.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>

#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

double msSince(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - t0)
      .count();
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "usage: bench_cage_h2h TRAJ [frame] [reps] [type]\n";
    return 2;
  }
  const std::string traj = argv[1];
  const int frame = argc > 2 ? std::atoi(argv[2]) : 1;
  const int reps = argc > 3 ? std::atoi(argv[3]) : 3;
  const int atomType = argc > 4 ? std::atoi(argv[4]) : 1;

  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO(traj, frame, yCloud, atomType);
  if (yCloud.nop == 0) {
    std::cerr << "could not read " << traj << " frame " << frame << "\n";
    return 1;
  }
  auto nList = nneigh::neighListO(3.5, yCloud, atomType);
  auto idx = nneigh::neighbourListByIndex(yCloud, nList);
  const auto ringsAll = primitive::ringNetwork(idx, 7);
  std::vector<std::vector<int>> six;
  for (const auto &r : ringsAll) {
    if (r.size() == 6) {
      six.push_back(r);
    }
  }
  const auto kept = ring::contractibleRings(ringsAll, yCloud);
  const auto rank = ring::positionRank(yCloud);
  std::vector<std::array<double, 3>> xyz(static_cast<size_t>(yCloud.nop));
  for (int i = 0; i < yCloud.nop; ++i) {
    xyz[static_cast<size_t>(i)] = {yCloud.pts[static_cast<size_t>(i)].x,
                                   yCloud.pts[static_cast<size_t>(i)].y,
                                   yCloud.pts[static_cast<size_t>(i)].z};
  }
  const auto sigs = std::vector<cage::Signature>{
      cage::Signature::parse("512"),
      cage::Signature::parse("51262"),
      cage::Signature::parse("51264"),
  };

  // one untimed pass for the counts
  std::vector<ring::strucType> ringType(six.size(), ring::strucType::unclassified);
  std::vector<cage::Cage> cages;
  const auto index = ring::buildRingSearchIndex(six, yCloud);
  auto listHC = ring::findHC(six, ringType, idx, cages, index);
  auto listDDC = ring::findDDC(six, ringType, listHC, cages, index);
  int nHC = 0;
  int nDDC = 0;
  int mixed = 0;
  int basal = 0;
  int prism = 0;
  ring::getStrucNumbers(ringType, cages, nHC, nDDC, mixed, prism, basal);
  int hcOnly = 0;
  int ddcOnly = 0;
  int both = 0;
  for (const auto t : ringType) {
    const bool h = t == ring::strucType::HCbasal || t == ring::strucType::HCprismatic ||
                   t == ring::strucType::bothBasal || t == ring::strucType::bothPrismatic;
    const bool d = t == ring::strucType::DDC || t == ring::strucType::bothBasal ||
                   t == ring::strucType::bothPrismatic;
    if (h && d) {
      ++both;
    } else if (h) {
      ++hcOnly;
    } else if (d) {
      ++ddcOnly;
    }
  }
  const auto found = cage::findBySignatures(kept, idx, sigs, rank, &xyz, &yCloud.box);
  const auto aff = ring::cageAffiliation(six, idx);

  const double cand = 3.5 + 1.5;
  const auto knn = nneigh::kNearestNeighbourPair(yCloud, 4, cand, atomType);
  const auto idxS = nneigh::neighbourListByIndex(yCloud, knn.first);
  const auto idxU = nneigh::neighbourListByIndex(yCloud, knn.second);
  auto sixOf = [](const std::vector<std::vector<int>> &all) {
    std::vector<std::vector<int>> sixOnly;
    for (const auto &r : all) {
      if (r.size() == 6) {
        sixOnly.push_back(r);
      }
    }
    return sixOnly;
  };
  const auto sixS = sixOf(primitive::ringNetwork(idxS, 6));
  const auto sixU = sixOf(primitive::ringNetwork(idxU, 6));
  const auto seeded =
      ring::seededCageAffiliation(sixS, idxS, sixU, idxU, false);
  int seedHc = 0;
  int seedDdc = 0;
  int seedBoth = 0;
  const int nAtom = std::min(yCloud.nop, static_cast<int>(seeded.hc.size()));
  for (int i = 0; i < nAtom; i++) {
    const bool h = seeded.hc[static_cast<size_t>(i)];
    const bool d = seeded.ddc[static_cast<size_t>(i)];
    if (h && d) {
      ++seedBoth;
    } else if (h) {
      ++seedHc;
    } else if (d) {
      ++seedDdc;
    }
  }

  auto fingerprint = [&](const std::vector<cage::Cage> &cl,
                         const std::vector<std::vector<cage::FoundCage>> &ff) {
    std::uint64_t hsh = 1469598103934665603ULL;
    auto mix = [&](std::uint64_t x) {
      hsh ^= x;
      hsh *= 1099511628211ULL;
    };
    for (const auto &c : cl) {
      mix(c.type == cage::cageType::HexC ? 1 : 2);
      for (const int r : c.rings) {
        mix(static_cast<std::uint64_t>(r + 3));
      }
      mix(0);
    }
    for (const auto &list : ff) {
      for (const auto &c : list) {
        for (const int f : c.faces) {
          mix(static_cast<std::uint64_t>(f + 5));
        }
        mix(1);
      }
      mix(2);
    }
    return hsh;
  };

  double bestIce = 1e300;
  double bestClath = 1e300;
  double bestBoth = 1e300;
  double bestAff = 1e300;
  double bestRings = 1e300;
  std::uint64_t finger = 0;
  for (int rep = 0; rep < reps; ++rep) {
    auto t0 = std::chrono::steady_clock::now();
    ringType.assign(six.size(), ring::strucType::unclassified);
    cages.clear();
    auto hc = ring::findHC(six, ringType, idx, cages, index);
    auto ddc = ring::findDDC(six, ringType, hc, cages, index);
    (void)ddc;
    bestIce = std::min(bestIce, msSince(t0));

    t0 = std::chrono::steady_clock::now();
    auto f = cage::findBySignatures(kept, idx, sigs, rank, &xyz, &yCloud.box);
    bestClath = std::min(bestClath, msSince(t0));

    t0 = std::chrono::steady_clock::now();
    ringType.assign(six.size(), ring::strucType::unclassified);
    cages.clear();
    hc = ring::findHC(six, ringType, idx, cages, index);
    ddc = ring::findDDC(six, ringType, hc, cages, index);
    f = cage::findBySignatures(kept, idx, sigs, rank, &xyz, &yCloud.box);
    const double bothMs = msSince(t0);
    bestBoth = std::min(bestBoth, bothMs);
    // One line per repetition so a median can be taken. sample_ms stays
    // the minimum, which is what earlier logs recorded.
    std::cout << "both_sample_ms " << bothMs << "\n";
    finger = fingerprint(cages, f);

    t0 = std::chrono::steady_clock::now();
    auto a = ring::cageAffiliation(six, idx);
    bestAff = std::min(bestAff, msSince(t0));
    (void)a;
    (void)f;

    t0 = std::chrono::steady_clock::now();
    auto rebuilt = primitive::ringNetwork(idx, 7);
    bestRings = std::min(bestRings, msSince(t0));
    (void)rebuilt;
  }

  std::cout << "atoms " << yCloud.nop << " rings " << ringsAll.size()
            << " contractible " << kept.size() << " six " << six.size() << "\n";
  std::cout << "cages HC " << nHC << " DDC " << nDDC << " mixedRings " << mixed
            << " basal " << basal << " prism " << prism << "\n";
  std::cout << "rings exclusive HC " << hcOnly << " DDC " << ddcOnly
            << " both " << both << " inclusive HC " << (hcOnly + both)
            << " DDC " << (ddcOnly + both) << "\n";
  std::cout << "listHC " << listHC.size() << " listDDC " << listDDC.size()
            << "\n";
  std::cout << "clathrate";
  for (size_t i = 0; i < found.size(); ++i) {
    std::cout << " " << sigs[i].str() << "=" << found[i].size();
  }
  int affHC = 0;
  int affDDC = 0;
  for (bool v : aff.hc) {
    affHC += v ? 1 : 0;
  }
  for (bool v : aff.ddc) {
    affDDC += v ? 1 : 0;
  }
  std::cout << "\naffiliation rings HC " << affHC << " DDC " << affDDC << "\n";
  std::cout << "seeded atoms HC-only " << seedHc << " DDC-only " << seedDdc
            << " both " << seedBoth << " HC-inclusive " << (seedHc + seedBoth)
            << " DDC-inclusive " << (seedDdc + seedBoth) << "\n";
  std::cout << "fingerprint " << std::hex << finger << std::dec << "\n";
  std::cout << "best_ms ice " << bestIce << " clathrate " << bestClath
            << " both " << bestBoth << " affiliation " << bestAff
            << " rings " << bestRings << "\n";
  std::cout << "sample_ms " << bestBoth << "\n";
  std::cout << "sample_rings_ms " << bestRings << "\n";
  return 0;
}
