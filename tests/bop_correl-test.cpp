// Internal
#include <bop_correl.hpp>

// Standard
#include <iostream>

// Conan
#include <catch2/catch.hpp>
// #include <rang.hpp>

// SCENARIO("Test the bond order correlation for BCC", "[bopMw]") {
//   GIVEN("An object for initializing bop") {
//     // DO NOT EDIT, FILE SPECIFIC
//     std::unique_ptr<chill::bop> bopTest(new chill::bop);
//     int atom_type = 1;
//     int frame_num = 1;
//     int num_of_neighbors = 4;
//     int num_of_points = 432;
//     std::array<double, 3> query_point = {0, 0, 0};
//     std::string testFile = "../input/traj/dump-bcc.lammpstrj";
//     chill::initSlice<double> starter;
//     chill::yodaPoint<double> testPoint;
//     // Fill starter
//     starter.coordHigh = {3, 3, 3};
//     starter.coordLow = {0, 0, 0};
//     starter.frameRange = {1, 2};
//     starter.filename = testFile;
//     WHEN("Initialze bop") {
//       bopTest->initBOP(num_of_points, atom_type, starter);
//       THEN("Something") {
//         // Hand calculate and test the Q_lm at pointQ(0)
//         // testPoint = bopTest->pointQ(0);
//         testPoint = bopTest->atomVerdict(0);
//         REQUIRE(1);
//       }
//     }
//     bopTest->cleanUp();
//   }
// }

SCENARIO("Test the bond order correlation for Mw Cubic", "[bopMw]") {
  GIVEN("An object for initializing bop") {
    // DO NOT EDIT, FILE SPECIFIC
    std::unique_ptr<chill::bop> bopTest(new chill::bop);
    int atom_type = 1;
    int frame_num = 1;
    int num_of_neighbors = 4;
    int num_of_points = 4096;
    std::array<double, 3> query_point = {0, 0, 0};
    std::string testFile = "../input/traj/mW_cubic.lammpstrj";
    chill::initSlice<double> starter;
    chill::yodaPoint<double> testPoint;
    chill::structurePercentage testPercent;
    // Fill starter
    // starter.coordHigh = {3, 3, 3};
    // starter.coordLow = {0, 0, 0};
    starter.frameRange = {1, 2};
    starter.filename = testFile;
    WHEN("Initialze bop") {
      bopTest->initBOP(num_of_points, atom_type, starter);
      THEN("Something") {
        // Hand calculate and test the Q_lm at pointQ(0)
        // testPoint = bopTest->pointQ(0);
        // testPoint = bopTest->atomVerdict(0);
        testPercent = bopTest->frameVerdict(1);
        REQUIRE(1);
      }
    }
    bopTest->cleanUp();
  }
}

// SCENARIO("Test the bond order correlation for Mw Hexagonal", "[bopMw]") {
//   GIVEN("An object for initializing bop") {
//     // DO NOT EDIT, FILE SPECIFIC
//     std::unique_ptr<chill::bop> bopTest(new chill::bop);
//     int atom_type = 1;
//     int frame_num = 1;
//     int num_of_neighbors = 4;
//     int num_of_points = 5096;
//     std::array<double, 3> query_point = {0, 2.597, 0.45};
//     std::string testFile = "../input/traj/mW_hexagonal.lammpstrj";
//     chill::initSlice<double> starter;
//     chill::yodaPoint<double> testPoint;
//     chill::structurePercentage testPercent;
//     // Fill starter
//     // starter.coordHigh = {3, 3, 3};
//     // starter.coordLow = {0, 0, 0};
//     starter.frameRange = {9, 10};
//     starter.filename = testFile;
//     WHEN("Initialze bop") {
//       bopTest->initBOP(num_of_points, atom_type, starter);
//       THEN("Something") {
//         // Hand calculate and test the Q_lm at pointQ(0)
//         // testPoint = bopTest->pointQ(0);
//         // testPercent = bopTest->frameVerdict(9);
//         REQUIRE(1);
//       }
//     }
//     bopTest->cleanUp();
//   }
// }
