//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

// Read a LAMMPS dump of a TIP4P/brine slab and print the in-situ front
// measurements. Ice labels are the published CHILL+ rule on oxygen type
// `--water-type`. Ions are `--ion-types`. Times are frame index times
// `--every` steps times `--dt-ps`.

#include <brine_front.hpp>
#include <chill_plus_scheme.hpp>
#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

bool crystal(molSys::atom_state_type label) {
  return label == molSys::atom_state_type::hexagonal || label == molSys::atom_state_type::cubic ||
         label == molSys::atom_state_type::reHex || label == molSys::atom_state_type::reCubic;
}

void sortById(Cloud &cloud) {
  std::sort(cloud.pts.begin(), cloud.pts.end(), gen::compareByAtomID);
  cloud.nop = static_cast<int>(cloud.pts.size());
  cloud.idIndexMap.clear();
  for (int i = 0; i < cloud.nop; ++i) {
    cloud.idIndexMap[cloud.pts[static_cast<std::size_t>(i)].atomID] = i;
  }
}

struct Args {
  std::string path;
  int waterType = 1;
  std::vector<int> ionTypes{3, 4};
  double cutoff = 3.5;
  int bins = 24;
  double dtPs = 2;
  int every = 5000;
  int timeRepeats = 0;
  int maxFrames = 0;
  std::string labels;
  bool profile = false;
};

std::vector<int> splitInts(const std::string &text) {
  std::vector<int> out;
  std::stringstream stream(text);
  std::string part;
  while (std::getline(stream, part, ',')) {
    if (!part.empty()) {
      out.push_back(std::stoi(part));
    }
  }
  return out;
}

Args parse(int argc, char **argv) {
  Args args;
  for (int i = 1; i < argc; ++i) {
    const std::string key = argv[i];
    auto need = [&](const char *name) {
      if (i + 1 >= argc) {
        throw std::runtime_error(std::string("missing value for ") + name);
      }
      return std::string(argv[++i]);
    };
    if (key == "--water-type") {
      args.waterType = std::stoi(need("--water-type"));
    } else if (key == "--ion-types") {
      args.ionTypes = splitInts(need("--ion-types"));
    } else if (key == "--cutoff") {
      args.cutoff = std::stod(need("--cutoff"));
    } else if (key == "--bins") {
      args.bins = std::stoi(need("--bins"));
    } else if (key == "--dt-ps") {
      args.dtPs = std::stod(need("--dt-ps"));
    } else if (key == "--every") {
      args.every = std::stoi(need("--every"));
    } else if (key == "--time") {
      args.timeRepeats = std::stoi(need("--time"));
    } else if (key == "--max-frames") {
      args.maxFrames = std::stoi(need("--max-frames"));
    } else if (key == "--labels") {
      args.labels = need("--labels");
    } else if (key == "--profile") {
      args.profile = true;
    } else if (key[0] != '-') {
      args.path = key;
    } else {
      throw std::runtime_error("unknown argument " + key);
    }
  }
  if (args.path.empty()) {
    throw std::runtime_error("usage: brine_front dump.lammpstrj [--time N]");
  }
  return args;
}

bool isIonType(const Args &args, int type) {
  return std::find(args.ionTypes.begin(), args.ionTypes.end(), type) != args.ionTypes.end();
}

struct Labelled {
  std::vector<bool> ice;
  std::vector<int> waters;
  std::vector<int> ions;
};

Labelled labelFrame(Cloud &cloud, const Args &args) {
  sortById(cloud);
  auto nList = nneigh::neighListO(args.cutoff, cloud, args.waterType);
  chill::applyChillPlus(cloud, nList, chill::ChillPlusRule::Published, false);
  Labelled labelled;
  labelled.ice.assign(static_cast<std::size_t>(cloud.nop), false);
  for (int i = 0; i < cloud.nop; ++i) {
    const auto &point = cloud.pts[static_cast<std::size_t>(i)];
    if (point.type == args.waterType) {
      labelled.waters.push_back(i);
      labelled.ice[static_cast<std::size_t>(i)] = crystal(point.iceType);
    } else if (isIonType(args, point.type)) {
      labelled.ions.push_back(i);
    }
  }
  return labelled;
}

double medianOf(std::vector<double> values) {
  if (values.empty()) {
    return 0;
  }
  std::sort(values.begin(), values.end());
  const std::size_t n = values.size();
  if (n % 2 == 1) {
    return values[n / 2];
  }
  return 0.5 * (values[n / 2 - 1] + values[n / 2]);
}

} // namespace

int main(int argc, char **argv) {
  try {
    const Args args = parse(argc, argv);
    site::FrontOptions options;
    options.nBins = args.bins;
    options.shellCutoff = args.cutoff;
    options.layerHalfWidth = args.cutoff;
    std::vector<double> timeNs;
    std::vector<double> thickness;
    std::vector<double> kValue;
    std::vector<double> cLiquid;
    int traps = 0;
    int releases = 0;
    site::IonEnvironment previous;
    bool havePrevious = false;
    int used = 0;
    Cloud timedCloud;
    bool haveTimed = false;
    for (int frame = 1; args.maxFrames <= 0 || used < args.maxFrames; ++frame) {
      Cloud cloud;
      sinp::readLammpsTrj(args.path, frame, cloud);
      if (cloud.nop <= 0 || cloud.pts.empty()) {
        break;
      }
      if (!haveTimed) {
        timedCloud = cloud;
        haveTimed = true;
      }
      auto labelled = labelFrame(cloud, args);
      const auto snap =
          site::analyseFront(cloud, labelled.ice, labelled.waters, labelled.ions, options);
      const double tNs = static_cast<double>(used) * args.every * args.dtPs / 1000.0;
      timeNs.push_back(tNs);
      thickness.push_back(snap.fronts.thickness);
      if (snap.rejection.defined) {
        kValue.push_back(snap.rejection.k);
      }
      cLiquid.push_back(snap.rejection.cLiquid);
      if (havePrevious) {
        for (const auto &event : site::trappingTransitions(previous, snap.ions)) {
          if (event.to == site::IonState::ice) {
            traps += 1;
          }
          if (event.from == site::IonState::ice) {
            releases += 1;
          }
        }
      }
      previous = snap.ions;
      havePrevious = true;
      if (!args.labels.empty() && used == 0) {
        std::ofstream labels(args.labels);
        labels << "# index atom x y z chill crystal\n";
        for (int water : labelled.waters) {
          const auto &point = cloud.pts[static_cast<std::size_t>(water)];
          labels << water << " " << point.atomID << " " << point.x << " " << point.y << " " << point.z
                 << " " << static_cast<int>(point.iceType) << " "
                 << (labelled.ice[static_cast<std::size_t>(water)] ? 1 : 0) << "\n";
        }
      }
      const int nCrystal = static_cast<int>(
          std::count(labelled.ice.begin(), labelled.ice.end(), true));
      std::cout << "frame " << used << " t_ns " << tNs << " n_crystal " << nCrystal
                << " thickness " << snap.fronts.thickness
                << " z_enter " << snap.fronts.zEnter << " z_exit " << snap.fronts.zExit
                << " c_ice " << snap.rejection.cIce << " c_liquid " << snap.rejection.cLiquid
                << " k " << snap.rejection.k << " k_defined " << (snap.rejection.defined ? 1 : 0)
                << " n_ice_ion " << snap.ions.nIce << " n_front_ion " << snap.ions.nFront
                << " n_liquid_ion " << snap.ions.nLiquid << " pockets " << snap.pockets.nPockets
                << " reservoir " << snap.pockets.reservoirWaters << " layer " << snap.layer.waters.size()
                << " defects " << snap.layer.nDefect << " four " << snap.layer.nFour
                << " rings6 " << (snap.layer.ringCensus.size() > 6 ? snap.layer.ringCensus[6] : 0)
                << " rings5 " << (snap.layer.ringCensus.size() > 5 ? snap.layer.ringCensus[5] : 0)
                << " hex " << snap.layer.nHexagonal << " cubic " << snap.layer.nCubic
                << " interfacial " << snap.layer.nInterfacial << "\n";
      if (!snap.pockets.sizes.empty()) {
        std::cout << "pocket_sizes";
        for (int size : snap.pockets.sizes) {
          std::cout << " " << size;
        }
        std::cout << "\n";
      }
      ++used;
    }
    const auto fit = site::linearFit(timeNs, thickness);
    const double front = fit.slope / 2.0;
    const double frontErr = fit.stderr / 2.0;
    double kMean = 0;
    for (double k : kValue) {
      kMean += k;
    }
    if (!kValue.empty()) {
      kMean /= static_cast<double>(kValue.size());
    }
    double cLast = cLiquid.empty() ? 0 : cLiquid.back();
    std::cout << "summary frames " << used << " front_A_per_ns " << front << " front_stderr " << frontErr
              << " thickness_slope " << fit.slope << " traps " << traps << " releases " << releases
              << " k_mean " << kMean << " c_liquid_last " << cLast << "\n";
    if (args.profile && haveTimed) {
      auto ms = [](auto a, auto b) {
        return std::chrono::duration<double, std::milli>(b - a).count();
      };
      const int repeats = args.timeRepeats >= 2 ? args.timeRepeats : 8;
      double neigh = 0, chill = 0, bins = 0, pockets = 0, rings = 0, ions = 0, total = 0;
      for (int i = 0; i < repeats; ++i) {
        Cloud again = timedCloud;
        auto t0 = std::chrono::steady_clock::now();
        std::sort(again.pts.begin(), again.pts.end(), gen::compareByAtomID);
        again.nop = static_cast<int>(again.pts.size());
        auto nList = nneigh::neighListO(args.cutoff, again, args.waterType);
        auto t1 = std::chrono::steady_clock::now();
        chill::applyChillPlus(again, nList, chill::ChillPlusRule::Published, false);
        auto t2 = std::chrono::steady_clock::now();
        std::vector<bool> ice(static_cast<std::size_t>(again.nop), false);
        std::vector<int> waters, ionIdx;
        for (int a = 0; a < again.nop; ++a) {
          const auto &point = again.pts[static_cast<std::size_t>(a)];
          if (point.type == args.waterType) {
            waters.push_back(a);
            ice[static_cast<std::size_t>(a)] = crystal(point.iceType);
          } else if (isIonType(args, point.type)) {
            ionIdx.push_back(a);
          }
        }
        auto prof = site::concentrationProfile(again, ice, waters, ionIdx, options.nBins, options.axis,
                                               options.ionsPerFormula);
        auto fronts = site::locateFronts(prof, options.frontCut);
        (void)site::rejectionCoefficient(prof, options.bulkIce, options.bulkLiquid);
        auto t3 = std::chrono::steady_clock::now();
        (void)site::brinePockets(again, ice, waters, ionIdx, options.shellCutoff);
        auto t4 = std::chrono::steady_clock::now();
        (void)site::frontLayer(again, waters, fronts, options.layerHalfWidth, options.shellCutoff,
                               options.maxRing, options.axis);
        auto t5 = std::chrono::steady_clock::now();
        (void)site::ionEnvironment(again, ice, ionIdx, 0, options.shellCutoff);
        auto t6 = std::chrono::steady_clock::now();
        neigh += ms(t0, t1);
        chill += ms(t1, t2);
        bins += ms(t2, t3);
        pockets += ms(t3, t4);
        rings += ms(t4, t5);
        ions += ms(t5, t6);
        total += ms(t0, t6);
      }
      const double n = repeats;
      std::cout << "profile_ms n " << repeats << " neigh " << neigh / n << " chill " << chill / n
                << " bins " << bins / n << " pockets " << pockets / n << " rings " << rings / n
                << " ions " << ions / n << " total " << total / n << "\n";
    }
    if (args.timeRepeats >= 2 && haveTimed) {
      std::vector<double> samples;
      samples.reserve(static_cast<std::size_t>(args.timeRepeats));
      for (int i = 0; i < args.timeRepeats; ++i) {
        Cloud again = timedCloud;
        const auto started = std::chrono::steady_clock::now();
        auto labelled = labelFrame(again, args);
        (void)site::analyseFront(again, labelled.ice, labelled.waters, labelled.ions, options);
        const auto finished = std::chrono::steady_clock::now();
        const double ms =
            std::chrono::duration<double, std::milli>(finished - started).count();
        samples.push_back(ms);
      }
      std::cout << "timing_ms n " << samples.size() << " median " << medianOf(samples);
      for (double sample : samples) {
        std::cout << " " << sample;
      }
      std::cout << "\n";
    }
    return used > 0 ? 0 : 1;
  } catch (const std::exception &error) {
    std::cerr << error.what() << "\n";
    return 1;
  }
}
