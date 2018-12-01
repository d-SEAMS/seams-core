///////////////////////////////////////////////////////////////////////////////////////////
//
//
// Harmonics Test
// MIT License

// Copyright (c) 2018 Amrita Goswami, Rohit Goswami
//
// amrita16thaug646[at]gmail.com, r95g10[at]gmail.com
// This code has been written for the purpose of obtaining the 1-D analogue of the Structure
// factor for a confined system, from a lammps trajectory file.
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
  GIVEN("A blaze static vector (2) of angles") {
    int orderL = 3;
    // {θ,ϕ}
    blaze::StaticVector<double, 2UL> testAngles{0, 1.047};
    blaze::StaticVector<std::complex<double>, 7UL> harmonicVector;
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
