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

#include <iostream>
#include <string>

// interface_order study <dump> [frame]
// interface_order growth <dump> <picoseconds between frames>
// interface_order q6 <dump> <picoseconds between frames>
// interface_order confined

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "usage: interface_order study <dump> [frame]\n"
                 "       interface_order growth <dump> <ps>\n"
                 "       interface_order q6 <dump> <ps>\n"
                 "       interface_order confined\n"
                 "       interface_order wc <dump> [frame]\n";
    return 2;
  }
  const std::string cmd = argv[1];
  if (cmd == "study" && argc >= 3) {
    const int frame = argc >= 4 ? std::stoi(argv[3]) : 1;
    const auto study = iface::studyDump(argv[2], frame);
    std::cout << iface::formatStudy(study, argv[2]);
    return study.n > 0 ? 0 : 1;
  }
  if (cmd == "growth" && argc >= 4) {
    const auto growth = iface::studyTrajectory(argv[2], std::stod(argv[3]));
    std::cout << iface::formatGrowth(growth);
    return growth.frames.empty() ? 1 : 0;
  }
  if (cmd == "q6" && argc >= 4) {
    const auto text = iface::formatQ6Growth(argv[2], std::stod(argv[3]));
    std::cout << text;
    return text.find("frame ") == std::string::npos ? 1 : 0;
  }
  if (cmd == "wc" && argc >= 3) {
    const int frame = argc >= 4 ? std::stoi(argv[3]) : 1;
    molSys::PointCloud<molSys::Point<double>, double> cloud;
    cloud = sinp::readLammpsTrjO(argv[2], frame, cloud, 1);
    std::cout << iface::formatWillardChandler(iface::willardChandler(std::move(cloud)));
    return 0;
  }
  if (cmd == "confined") {
    std::cout << iface::formatStudy(iface::studyCloud(iface::nanotube(4, 8)), "square");
    std::cout << iface::formatStudy(iface::studyCloud(iface::nanotube(5, 8)), "pentagon");
    std::cout << iface::formatStudy(iface::studyCloud(iface::nanotube(6, 8)), "hexagon");
    std::cout << iface::formatStudy(iface::studyCloud(iface::honeycombSheet(4, 3, 1)),
                                    "monolayer");
    std::cout << iface::formatStudy(iface::studyCloud(iface::honeycombSheet(4, 3, 2)),
                                    "bilayer");
    return 0;
  }
  std::cerr << "unknown command\n";
  return 2;
}
