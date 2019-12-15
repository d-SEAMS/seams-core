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
    int nop = templatePntCloud.nop;             // Number of particles in an HC
    int dim = 3;                                // Number of dimensions
    // Eigen matrix for the reference (or right pointSet)
    Eigen::MatrixXd refPointSet(nop, dim);
    // Read in the XYZ template file
    sinp::readXYZ(refFileXYZ, &templatePntCloud);
    // Read in the point-to-point correspondence vector of vectors
    refPntToPnt = sinp::readRefHCdata(refFileData);
    // Get the reference point set
    refPointSet = pntToPnt::fillPointSetHC(&templatePntCloud, refPntToPnt, 0);
    // --------------------------
    // GETTING THE ROTATED POINT SET
    std::vector<std::vector<int>>
        targetPntToPnt;  // Vector of vector of ints with
                         // the connectivity information for the target point
                         // set
    molSys::PointCloud<molSys::Point<double>, double>
        targetCloud;  // pointCloud
    molSys::Point<double> iPoint;
    Eigen::MatrixXd targetPointSet(nop, dim);
    // Fill the pointCloud
    //
    iPoint.type = 1;  // Same for all the points here
    // Element {0}
    iPoint.atomID = 0;         // iatom
    iPoint.x = 19.8443861192;  // x
    iPoint.y = 3.6525394722;   // y
    iPoint.z = 15.0939999;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {1}
    iPoint.atomID = 1;         // iatom
    iPoint.x = 15.9491946129;  // x
    iPoint.y = 5.9012086664;   // y
    iPoint.z = 15.0939999;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {2}
    iPoint.atomID = 2;         // iatom
    iPoint.x = 15.9490040262;  // x
    iPoint.y = 1.4035395722;   // y
    iPoint.z = 15.0939999;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {3}
    iPoint.atomID = 3;         // iatom
    iPoint.x = 18.5458860692;  // x
    iPoint.y = 5.9016075325;   // y
    iPoint.z = 14.1949997;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {4}
    iPoint.atomID = 4;         // iatom
    iPoint.x = 18.5456945129;  // x
    iPoint.y = 1.4039389177;   // y
    iPoint.z = 14.1949997;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {5}
    iPoint.atomID = 5;         // iatom
    iPoint.x = 14.6505039762;  // x
    iPoint.y = 3.6526076325;   // y
    iPoint.z = 14.1949997;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {6}
    iPoint.atomID = 6;         // iatom
    iPoint.x = 18.5458860692;  // x
    iPoint.y = 5.9016075325;   // y
    iPoint.z = 11.4329996;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {7}
    iPoint.atomID = 7;         // iatom
    iPoint.x = 18.5456945129;  // x
    iPoint.y = 1.4039389177;   // y
    iPoint.z = 11.4329996;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {8}
    iPoint.atomID = 8;         // iatom
    iPoint.x = 14.6505039762;  // x
    iPoint.y = 3.6526076325;   // y
    iPoint.z = 11.4329996;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {9}
    iPoint.atomID = 9;         // iatom
    iPoint.x = 19.8443861192;  // x
    iPoint.y = 3.6525394722;   // y
    iPoint.z = 10.5340004;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {10}
    iPoint.atomID = 10;        // iatom
    iPoint.x = 15.9491946129;  // x
    iPoint.y = 5.9012086664;   // y
    iPoint.z = 10.5340004;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {11}
    iPoint.atomID = 11;        // iatom
    iPoint.x = 15.9490040262;  // x
    iPoint.y = 1.4035395722;   // y
    iPoint.z = 10.5340004;     // z
    targetCloud.pts.push_back(iPoint);
    // Update nop
    targetCloud.nop = targetCloud.pts.size();
    //
    //
    targetPntToPnt = refPntToPnt;  // In this case, they are the same
    //
    //
    // Get the target point set
    targetPointSet = pntToPnt::fillPointSetHC(&targetCloud, targetPntToPnt, 0);
    // --------------------------
    // Now get the absolute orientation of the left (candidate/target) system
    // with respect to the right (template/reference) system test
    //
    std::vector<double> quaternionRot;  // quaternion rotation
    //
    // absolute orientation using Horn's algorithm between the target and test
    // set
    absor::hornAbsOrientation(refPointSet, targetPointSet, &quaternionRot);
    //
    std::vector<double>
        selfQuatRot;  // quaternion for the reference set and itself
    absor::hornAbsOrientation(refPointSet, refPointSet, &selfQuatRot);

    double angDist = gen::angDistDegQuaternions(selfQuatRot, quaternionRot);
    //
    REQUIRE_THAT(angDist, Catch::Matchers::Floating::WithinAbsMatcher(
                              30.0, 0.1));  // Evaluate condition
    // --------------------------
  }  // End of given
}  // End of scenario