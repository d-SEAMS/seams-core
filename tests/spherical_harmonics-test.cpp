///////////////////////////////////////////////////////////////////////////////////////////
//
//
// Harmonics Test
// d-SEAMS molecular dynamics analysis engine code
// Copyright (C) <2018-present> Amrita Goswami, Rohit Goswami
// amrita16thaug646[at]gmail.com, r95g10[at]gmail.com

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
//
///////////////////////////////////////////////////////////////////////////////////////////
// Internal
#include "spherical_harmonics.h"

// Standard
#include <array>
#include <iostream>

// Conan
#include <catch2/catch.hpp>
#include <rang.hpp>

SCENARIO("Get θ and ϕ from cartesian coordinates", "[sphereAngle]") {
  GIVEN("A cartesian coordinate array.") {
    std::array<double, 3> testPoint{1.732, 0, 1};
    std::array<double, 2> convertedPoint{0, 0};
    WHEN("A coordinate transform occurs") {
      convertedPoint = trans::radialCoord(testPoint);
      THEN("We get the polar and azimuthal angles.") {
        REQUIRE(convertedPoint.size() == 2);
        REQUIRE_THAT(convertedPoint[0], Catch::Matchers::WithinAbs(0, 1.0e-10));
        REQUIRE_THAT(convertedPoint[1],
                     Catch::Matchers::WithinAbs(1.047184849, 1.0e-10));
        // Test
        std::cout << std::endl;
        std::cout << rang::style::bold << "<x,y,z>" << std::endl
                  << rang::style::reset;
        for (const auto &s : testPoint)
          std::cout << rang::fg::blue << s << ' ' << rang::style::reset;

        std::cout << std::endl;
        std::cout << rang::style::bold << "<ϕ,θ>" << std::endl
                  << rang::style::reset;
        for (const auto &s : convertedPoint)
          std::cout << rang::fg::green << s << ' ' << rang::style::reset;
        std::cout << std::endl;
      }
    }
  }
}

// TODO: Test transposed result vectors by value, like in the neighbors test
SCENARIO("Test the boost spherical harmonics", "[sphericalHarmonics]") {
  GIVEN("A static vector (2) of angles") {
    int orderL = 3;
    // {θ,ϕ}
    std::array<double, 2> testAngles{0, 1.047};
    std::vector<std::complex<double>> harmonicVector;
    WHEN("Yₗₘ is calculated") {
      harmonicVector = trans::spheriHarmo(orderL, testAngles);
      THEN("We get the polar and azimuthal angles.") {
        REQUIRE(harmonicVector.size() == 7);
        // Test
        std::cout << std::endl;
        std::cout << rang::style::bold << "<θ,ϕ>" << std::endl
                  << rang::style::reset;
        for (const auto &s : testAngles)
          std::cout << rang::fg::green << s << ' ' << rang::style::reset;
        std::cout << std::endl;
        std::cout << rang::style::bold << "Transposed result vector"
                  << std::endl
                  << rang::style::reset;
        for (const auto &s : harmonicVector)
          std::cout << rang::fg::blue << s << ' ' << rang::style::reset;
        std::cout << std::endl;
      }
    }
  }
}

TEST_CASE("No negative values for harmonics", "[radialNegative]") {
  std::srand(std::time(0)); //use current time as seed for random generator
  // Initialize Variables
  double random_variable = std::rand();
  std::array<double, 3> testPoint;
  std::array<double, 2> null{0, 0};
  // Always Positive
  testPoint = {random_variable, random_variable, random_variable};
  REQUIRE(trans::radialCoord(testPoint) >= null);
  testPoint = {-random_variable, random_variable, random_variable};
  REQUIRE(trans::radialCoord(testPoint) >= null);
  testPoint = {-random_variable, random_variable, -random_variable};
  REQUIRE(trans::radialCoord(testPoint) >= null);
  // Always Negative
  testPoint = {-random_variable, -random_variable, random_variable};
  REQUIRE(trans::radialCoord(testPoint) <= null);
  testPoint = {-random_variable, -random_variable, -random_variable};
  REQUIRE(trans::radialCoord(testPoint) <= null);
}

TEST_CASE("Test wraps", "[testwrap]") {
  // Initialize Variables
  std::array<double, 3> testPoint{1.592, 1.592, 1.592};
  std::array<double, 3> testPoint1{-1.592, -1.592, -1.592};
  std::array<double, 2> null{0, 0};
  REQUIRE(trans::radialCoord(testPoint) > null);
  REQUIRE(trans::radialCoord(testPoint1) < null);
  std::array<double, 2> testAngles = trans::radialCoord(testPoint);
  std::vector<std::complex<double>> harmonicVector =
      trans::spheriHarmo(3, testAngles);
  REQUIRE(trans::radialCoord(testPoint1) < null);
  std::array<double, 2> testAngles1 = trans::radialCoord(testPoint1);
  std::vector<std::complex<double>> harmonicVector1 =
      trans::spheriHarmo(3, testAngles1);
  REQUIRE(1);
}

TEST_CASE("Integration test for Ylm", "[intylm]") {
  // Initialize Variables
  std::array<double, 3> testPoint{0, -3.185, -3.185};
  std::array<double, 2> null{0, 0};
  REQUIRE(trans::radialCoord(testPoint) < null);
  std::array<double, 2> testAngles = trans::radialCoord(testPoint);
  std::vector<std::complex<double>> harmonicVector =
      trans::spheriHarmo(3, testAngles);
  REQUIRE(1);
}

TEST_CASE("Integration test for hexagonal mW", "[inthex]") {
  // Initialize Variables
  std::array<double, 3> testPoint{-2.249, 1.299, -3.66};
  std::array<double, 2> null{0, 0};
  REQUIRE(trans::radialCoord(testPoint) > null);
  std::array<double, 2> testAngles = trans::radialCoord(testPoint);
  std::vector<std::complex<double>> harmonicVector =
      trans::spheriHarmo(3, testAngles);
  REQUIRE(1);
}
