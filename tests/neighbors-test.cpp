// Internal
#include <neighbors.hpp>

// Standard
#include <iostream>

// Conan
#include <catch2/catch.hpp>
#include <rang.hpp>

SCENARIO("Test the neighborlist (number) generation", "[KNNneighborlist]") {
  GIVEN("A file with coordinates and the number of nearest neighbors") {
    // DO NOT EDIT, FILE SPECIFIC
    auto tkn = new neigh::treeKNN;
    std::array<double, 3> testBox = {};
    neigh::PointCloud<double> resultCloud;
    int atom_type = 1;
    int frame_num = 1;
    int num_of_neighbors = 5;
    int num_of_points = 432;
    std::string filename = "../input/traj/dump-bcc.lammpstrj";
    WHEN("Cloud is populated for a type") {
      tkn->initKNN(num_of_points, filename, frame_num, atom_type);
      THEN("We run the KNN search") {
        resultCloud = tkn->byNumber(0, num_of_neighbors);
        THEN("We get a resultCloud") {
          THEN("Which should the same size as the number of neighbors") {
            REQUIRE(resultCloud.pts.size() == num_of_neighbors);
          }
          THEN("We get coordinates of the nearest neighbors") {
            // Nearest neighbor
            REQUIRE(resultCloud.pts[0].x == 0.62);
            REQUIRE(resultCloud.pts[0].y == 0.62);
            REQUIRE(resultCloud.pts[0].z == 0.62);
            // Second neighbor
            REQUIRE(resultCloud.pts[1].x == 0);
            REQUIRE(resultCloud.pts[1].y == 0);
            REQUIRE(resultCloud.pts[1].z == 1.24);
            // Third neighbor
            REQUIRE(resultCloud.pts[2].x == 1.24);
            REQUIRE(resultCloud.pts[2].y == 0);
            REQUIRE(resultCloud.pts[2].z == 0);
            // Fourth neighbor
            REQUIRE(resultCloud.pts[3].x == 0);
            REQUIRE(resultCloud.pts[3].y == 1.24);
            REQUIRE(resultCloud.pts[3].z == 0);
            // Fifth neighbor
            REQUIRE(resultCloud.pts[4].x == 1.24);
            REQUIRE(resultCloud.pts[4].y == 1.24);
            REQUIRE(resultCloud.pts[4].z == 0);
          }
          THEN("We also get the box dimensions") {
            REQUIRE_THAT(resultCloud.box[0], Catch::Matchers::WithinAbs(
                                                 7.4379999999999997, 1.0e-10));
            REQUIRE_THAT(resultCloud.box[1], Catch::Matchers::WithinAbs(
                                                 7.4379999999999997, 1.0e-10));
            REQUIRE_THAT(resultCloud.box[2], Catch::Matchers::WithinAbs(
                                                 7.4379999999999997, 1.0e-10));
          }
        }
      }
    }
  }
}

// TODO: Pretty print like the sph harmonics test
