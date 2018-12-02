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
    int atom_type = 1;
    int frame_num = 1;
    int num_of_neighbors = 5;
    int num_of_points = 432;
    std::string testFile = "../input/traj/dump-bcc.lammpstrj";
    chill::initSlice<double> starter;
    // Fill starter
    starter.coordHigh = {3, 3, 3};
    starter.coordLow = {1, 1, 1};
    starter.frameRange = {1, 2};
    starter.filename = testFile;
    WHEN("Initialze bop") {
      bopTest->initBOP(starter);
      THEN("SOMETHING") { REQUIRE(1); }
    }
  }
}
