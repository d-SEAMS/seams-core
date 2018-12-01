// Internal
#include <neighbors.hpp>

// Standard
#include <iostream>

// Conan
#include <catch2/catch.hpp>
#include <rang.hpp>

SCENARIO("Test the neighborlist (number) generation", "[KNNneighborlist]") {
  GIVEN("A file with coordinates and the number of nearest neighbors") {
    auto tkn = new neigh::treeKNN;
    std::array<double, 3> testBox;
    // std::array<double, 3> cH{3.2, 5.2, 5.2};
    // std::array<double, 3> cL{0, 0, 0};
    neigh::PointCloud<double> resultCloud;
    int atom_type = 1;
    int frame_num = 1;
    int num_of_points = 432;
    std::string filename = "../input/traj/dump-bcc.lammpstrj";
    WHEN("Cloud is populated for a type") {
      tkn->initKNN(num_of_points, filename, frame_num, atom_type);
      THEN("We run the KNN search") {
        resultCloud = tkn->byNumber(0, 5);
        THEN("We get a resultCloud") {
          REQUIRE(1);
          // // Test
          // std::cout << std::endl;
          // std::cout << rang::style::bold << "<θ,ϕ>" << std::endl
          //           << rang::style::reset;
          // for (const auto &s : testAngles)
          //   std::cout << rang::fg::green << s << ' ' << rang::style::reset;
          // std::cout << std::endl;
          // std::cout << rang::style::bold << "Transposed result vector"
          //           << std::endl
          //           << rang::style::reset;
          // for (const auto &s : harmonicVector)
          //   std::cout << rang::fg::blue << s << ' ' << rang::style::reset;
          // std::cout << std::endl;
        }
      }
    }
  }
}

// TEST_CASE("Test nanoflann", "[nanoflann]") {
//   // srand(static_cast<unsigned int>(time(nullptr)));
//   // kdtree_demo<float>(1000000);
//   auto tkn = new neigh::treeKNN;
//   REQUIRE(1 >= 0);
// }
