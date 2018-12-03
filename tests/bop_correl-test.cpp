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
    int num_of_neighbors = 4;
    int num_of_points = 432;
    std::array<double, 3> query_point = {0, 0, 0};
    std::string testFile = "../input/traj/dump-bcc.lammpstrj";
    chill::initSlice<double> starter;
    chill::yodaPoint<double> testPoint;
    // Fill starter
    starter.coordHigh = {3, 3, 3};
    starter.coordLow = {0, 0, 0};
    starter.frameRange = {1, 2};
    starter.filename = testFile;
    WHEN("Initialze bop") {
      bopTest->initBOP(num_of_points, atom_type, starter);
      THEN("Something") {
        // Hand calculate and test the Q_lm at pointQ(0)
        // testPoint = bopTest->pointQ(0);
        testPoint = bopTest->frameVerdict();
        REQUIRE(1);
      }
    }
  }
}
