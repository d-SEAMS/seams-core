//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//-----------------------------------------------------------------------------------

#include <interface_order.hpp>
#include <seams_input.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using Catch::Approx;

TEST_CASE("a Gaussian-smoothed step has the normal-cdf width", "[iface]") {
  REQUIRE(iface::gaussianStepWidth(2.4) == Approx(2.0 * 1.2815515655446004 * 2.4));
  REQUIRE(iface::gaussianStepWidth(1.0) == Approx(2.5631031310892008));
}

TEST_CASE("bootstrap intervals on a paired difference", "[iface]") {
  const std::vector<double> sameA = {1.0, 2.0, 3.0, 4.0};
  const std::vector<double> sameB = {1.0, 2.0, 3.0, 4.0};
  const auto none = iface::pairedDifference(sameA, sameB, 2000, 270);
  REQUIRE(none.mean == Approx(0.0));
  REQUIRE(none.excludesZero == false);
  const std::vector<double> high = {10.0, 10.0, 10.0, 10.0};
  const std::vector<double> low = {0.0, 0.0, 0.0, 0.0};
  const auto gap = iface::pairedDifference(high, low, 500, 270);
  REQUIRE(gap.mean == Approx(10.0));
  REQUIRE(gap.excludesZero);
  REQUIRE(gap.lo > 0.0);
}

TEST_CASE("a straight line has slope two", "[iface]") {
  const std::vector<double> x = {0.0, 1.0, 2.0, 3.0};
  const std::vector<double> y = {10.0, 12.0, 14.0, 16.0};
  REQUIRE(iface::slope(x, y) == Approx(2.0));
  REQUIRE(iface::slope({}, {}) == 0.0);
}

TEST_CASE("ideal ice Ih is the same crystal under both CHILL+ rules", "[iface]") {
  const auto ice = iface::studyDump("traj/ice1h.lammpstrj", 1);
  REQUIRE(ice.n == 128);
  REQUIRE(ice.nDisagree == 0);
  REQUIRE(ice.countPub[1] == 128);
  REQUIRE(ice.countOvi[1] == 128);
  REQUIRE(ice.cnaCrystal == 0);
}

TEST_CASE("the 2015 rule and the OVITO rule disagree on the published counts",
          "[iface]") {
  const auto nucleus = iface::studyDump("traj/nucleation_f1.lammpstrj", 1);
  INFO(iface::formatStudy(nucleus, "nucleation"));
  REQUIRE(nucleus.n == 4096);
  REQUIRE(nucleus.nDisagree == 397);
  int causes = 0;
  for (int i = 1; i <= 4; ++i) {
    causes += nucleus.causeCount[i];
  }
  REQUIRE(causes == 397);
  REQUIRE(nucleus.causeCount[0] == 4096 - 397);

  const auto face = iface::studyDump("traj/interface_f0.lammpstrj", 1);
  INFO(iface::formatStudy(face, "interface"));
  REQUIRE(face.n == 864);
  REQUIRE(face.nDisagree == 46);
  causes = 0;
  for (int i = 1; i <= 4; ++i) {
    causes += face.causeCount[i];
  }
  REQUIRE(causes == 46);
  REQUIRE(face.kind == iface::Kind::slab);
  REQUIRE(face.solidPub.defined);
  REQUIRE(face.solidOvi.defined);
  REQUIRE(face.depthDisagree < face.depthAgree);

  molSys::PointCloud<molSys::Point<double>, double> nucleusCloud;
  nucleusCloud = sinp::readLammpsTrjO("traj/nucleation_f1.lammpstrj", 1, nucleusCloud, 1);
  const auto surface = iface::willardChandler(std::move(nucleusCloud));
  INFO(iface::formatWillardChandler(surface));
  REQUIRE(surface.defined);
  REQUIRE(surface.densityWidth > 0.0);
  REQUIRE(surface.densityWidth < nucleus.solidPub.fronts.at(0).width);
  REQUIRE(surface.densityWidth < nucleus.solidOvi.fronts.at(0).width);
}

TEST_CASE("prism blocks and sheet rings name confined ice that CHILL+ CNA and PTM do not",
          "[iface]") {
  const auto square = iface::studyCloud(iface::nanotube(4, 8));
  INFO(iface::formatStudy(square, "square"));
  REQUIRE(square.n == 32);
  REQUIRE(square.prismBlocks > 0);
  REQUIRE(square.prismAtoms == square.n);
  REQUIRE(square.bilayerAtoms == 0);
  REQUIRE(square.cnaCrystal == 0);
  REQUIRE(square.diamondCubic + square.diamondHex == 0);
  REQUIRE(square.countPub[0] + square.countPub[1] < square.n);

  const auto pentagon = iface::studyCloud(iface::nanotube(5, 8));
  INFO(iface::formatStudy(pentagon, "pentagon"));
  REQUIRE(pentagon.prismBlocks > 0);
  REQUIRE(pentagon.cnaCrystal == 0);
  REQUIRE(pentagon.ptmPass == 0);

  const auto hex = iface::studyCloud(iface::nanotube(6, 8));
  INFO(iface::formatStudy(hex, "hexagon"));
  REQUIRE(hex.prismBlocks > 0);
  REQUIRE(hex.bilayerAtoms == 0);
  REQUIRE(hex.cnaCrystal == 0);

  const auto mono = iface::studyCloud(iface::honeycombSheet(4, 3, 1));
  INFO(iface::formatStudy(mono, "monolayer"));
  REQUIRE(mono.thinSheets == 1);
  REQUIRE(mono.polygonAtoms > 0);
  REQUIRE(mono.prismBlocks == 0);
  REQUIRE(mono.bilayerAtoms == 0);
  REQUIRE(mono.countPub[2] == mono.n);
  REQUIRE(mono.cnaCrystal == 0);
  REQUIRE(mono.ptmPass == 0);

  const auto bilayer = iface::studyCloud(iface::honeycombSheet(4, 3, 2));
  INFO(iface::formatStudy(bilayer, "bilayer"));
  REQUIRE(bilayer.thinSheets == 2);
  REQUIRE(bilayer.bilayerAtoms > 0);
  REQUIRE(bilayer.prismBlocks == 0);
  REQUIRE(bilayer.cnaCrystal == 0);
  REQUIRE(bilayer.diamondCubic + bilayer.diamondHex == 0);
  REQUIRE(bilayer.countPub[0] + bilayer.countPub[1] < bilayer.n);
}
