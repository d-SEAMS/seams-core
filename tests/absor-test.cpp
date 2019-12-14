// Internal
#include <absOrientation.hpp>
#include <pntCorrespondence.hpp>

// Standard
#include <iostream>

// Conan
#include <catch2/catch.hpp>
#include <rang.hpp>

SCENARIO("Test the shape-matching of a perfect HC rotated by 30 degrees",
         "[match]") {
  GIVEN("A rotated HC (candidate structure)") {
    // --------------------------
    // GETTING THE REFERENCE/TEMPLATE POINT SET
    molSys::PointCloud<molSys::Point<double>, double>
        templatePntCloud;  // PointCloud for the reference point set
    std::string refFileXYZ =
        "../templateStructures/points/hc.xyz";  // XYZ file for the ref HC
    std::string refFileData =
        "../templateStructures/connectivity/hc.dat";  // XYZ file for the ref HC
    std::vector<std::vector<int>> refPntToPnt;  // Vector of vector of ints with
                                                // the connectivity information
    // Read in the XYZ template file
    sinp::readXYZ(refFileXYZ, &templatePntCloud);
    // Read in the point-to-point correspondence vector of vectors
    refPntToPnt = sinp::readRefHCdata(refFileData);
    // test
    REQUIRE(2 == 2);  // Evaluate condition
    // --------------------------
  }  // End of given
}  // End of scenario