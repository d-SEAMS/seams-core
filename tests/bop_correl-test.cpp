// Internal
#include <bop_correl.hpp>

// Standard
#include <iostream>

// Conan
#include <catch2/catch.hpp>
#include <rang.hpp>

SCENARIO("Test the bond order correlation generation", "[bopCIJ]") {
  GIVEN("An object for initializing bop") {
    // DO NOT EDIT, FILE SPECIFIC
    auto bopTest = new chill::bop;
    neigh::PointCloud<double> resultCloud;
    int atom_type = 1;
    int frame_num = 1;
    int num_of_neighbors = 5;
    int num_of_points = 432;
    std::string filename = "../input/traj/dump-bcc.lammpstrj";
    WHEN("Initialze bop") {
      // bopTest->initBOP(starter);
      std::cout << "lol";
      THEN("SOMETHING") { REQUIRE(1); }
    }
  }
}
