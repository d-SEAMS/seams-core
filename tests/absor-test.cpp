// Internal
#include <absOrientation.hpp>
#include <bulkTUM.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <pntCorrespondence.hpp>
#include <ring.hpp>
#include <topo_bulk.hpp>

// Standard
#include <iostream>

#include <catch2/catch.hpp>
#include <rang.hpp>

SCENARIO("Test the shape-matching of a perfect HC rotated by 30 degrees",
         "[match]") {
  GIVEN("A rotated HC (candidate structure)") {
    // --------------------------
    // GETTING THE REFERENCE/TEMPLATE POINT SET
    int nop = 12;  // Number of particles in an HC
    int dim = 3;   // Number of dimensions
    //
    std::string filePathXYZ = "../templates/hc.xyz";
    // Variables for rings
    std::vector<std::vector<int>> nList;  // Neighbour list
    std::vector<std::vector<int>> rings;  // Rings
    std::vector<ring::strucType>
        ringType;  // This vector will have a value for each ring inside
    std::vector<int> listHC;  // Contains atom indices of atoms making up HCs
    // Make a list of all the DDCs and HCs
    std::vector<cage::Cage> cageList;
    Eigen::MatrixXd refPnts(12, 3);  // Reference point set (Eigen matrix)
    int iring, jring;
    // -------------------------------------
    //
    // REFERENCE POINT SET
    refPnts = tum3::buildRefHC(filePathXYZ);
    //
    // -------------------------------------
    // GETTING THE ROTATED POINT SET
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
    // box lengths
    targetCloud.box.push_back(50);  // x box length
    targetCloud.box.push_back(50);  // y box length
    targetCloud.box.push_back(50);  // z box length
    // Update the unordered map
    for (int iatom = 0; iatom < targetCloud.nop; iatom++) {
      targetCloud.idIndexMap[iatom] = iatom;
    }  // end of filling the map
    //
    //
    // --------------------------
    // GETTING THE TARGET POINT SET
    // Calculate a neighbour list
    nList = nneigh::neighListO(3.5, &targetCloud, 1);
    // Neighbour list by index
    nList = nneigh::neighbourListByIndex(&targetCloud, nList);
    // Find the vector of vector of rings
    rings = primitive::ringNetwork(nList, 6);
    // init the ringType vector
    ringType.resize(rings.size());
    // Find the HCs
    listHC = ring::findHC(rings, &ringType, nList, &cageList);
    // Get the basal rings from cageList
    iring = cageList[0].rings[0];
    jring = cageList[0].rings[1];
    //
    std::vector<int> matchedBasal1,
        matchedBasal2;  // Re-ordered basal rings 1 and 2
    // Reordered basal rings
    // Getting the target Eigen vectors
    // Get the re-ordered matched basal rings, ordered with respect to each
    // other

    pntToPnt::relOrderHC(&targetCloud, rings[iring], rings[jring], nList,
                         &matchedBasal1, &matchedBasal2);
    //
    // --------------------------
    // Now get the absolute orientation of the left (candidate/target) system
    // with respect to the right (template/reference) system test
    //
    std::vector<double> quaternionRot;         // quaternion rotation
    double rmsd1, rmsd2;                       // least RMSD
    std::vector<double> rmsdList1, rmsdList2;  // List of RMSD per atom
    double scale;                              // Scale factor
    //
    //
    // Variables for looping through possible permutations
    //
    std::vector<double> currentQuat;      // quaternion rotation
    double currentRmsd;                   // least RMSD
    std::vector<double> currentRmsdList;  // List of RMSD per atom
    double currentScale;
    // absolute orientation using Horn's algorithm between the target and test
    // set
    int index;
    for (int i = 0; i < 6; i++) {
      // Change the order of the target points somehow!
      //
      targetPointSet = pntToPnt::changeHexCageOrder(&targetCloud, matchedBasal1,
                                                    matchedBasal2, i);
      // Shape-matching
      absor::hornAbsOrientation(refPnts, targetPointSet, &currentQuat,
                                &currentRmsd, &currentRmsdList, &currentScale);
      if (i == 0) {
        quaternionRot = currentQuat;
        rmsd1 = currentRmsd;
        rmsdList1 = currentRmsdList;
        scale = currentScale;
        index = 0;
      } else {
        if (currentRmsd < rmsd1) {
          quaternionRot = currentQuat;
          rmsd1 = currentRmsd;
          rmsdList1 = currentRmsdList;
          scale = currentScale;
          index = i;
        }  // update
      }    // Update if this is a better match
    }      // Loop through possible permutations
    // ---------
    std::vector<double>
        selfQuatRot;   // quaternion for the reference set and itself
    double selfScale;  // Scale for the reference set and itself

    // Shape-matching
    absor::hornAbsOrientation(refPnts, refPnts, &selfQuatRot, &rmsd2,
                              &rmsdList2, &selfScale);

    //
    double angDist = gen::angDistDegQuaternions(selfQuatRot, quaternionRot);
    //
    REQUIRE_THAT(angDist, Catch::Matchers::Floating::WithinAbsMatcher(
                              30.0, 0.01));  // Evaluate condition
    // --------------------------
  }  // End of given
}  // End of scenario

// DDC shape-matching
SCENARIO("Test the shape-matching of a perfect DDC rotated by 30 degrees",
         "[match]") {
  GIVEN("A rotated DDC (candidate structure)") {
    // --------------------------
    // GETTING THE REFERENCE/TEMPLATE POINT SET
    int nop = 14;  // Number of particles in an HC
    int dim = 3;   // Number of dimensions
    //
    std::string filePathXYZ = "../templates/ddc.xyz";
    // Variables for rings
    std::vector<std::vector<int>> nList;  // Neighbour list
    std::vector<std::vector<int>> rings;  // Rings
    std::vector<ring::strucType>
        ringType;  // This vector will have a value for each ring inside
    std::vector<int> listHC,
        listDDC;  // Contains atom indices of atoms making up HCs
    // Make a list of all the DDCs and HCs
    std::vector<cage::Cage> cageList;
    std::vector<int> ddcOrder;       // Connectivity of the DDC
    Eigen::MatrixXd refPnts(14, 3);  // Reference point set (Eigen matrix)
    // -------------------------------------
    //
    // REFERENCE POINT SET
    refPnts = tum3::buildRefDDC(filePathXYZ);
    //
    // -------------------------------------
    // GETTING THE ROTATED POINTS
    molSys::PointCloud<molSys::Point<double>, double>
        targetCloud;  // pointCloud
    molSys::Point<double> iPoint;
    Eigen::MatrixXd targetPointSet(nop, dim);
    // Fill the pointCloud
    //
    iPoint.type = 1;  // Same for all the points here
    // Element {0}
    iPoint.atomID = 0;      // iatom
    iPoint.x = 9.3263264;   // x
    iPoint.y = 34.8063278;  // y
    iPoint.z = 19.1100006;  // z
    targetCloud.pts.push_back(iPoint);
    // Element {1}
    iPoint.atomID = 1;      // iatom
    iPoint.x = 9.3263264;   // x
    iPoint.y = 34.8063278;  // y
    iPoint.z = 25.4799996;  // z
    targetCloud.pts.push_back(iPoint);
    // Element {2}
    iPoint.atomID = 2;      // iatom
    iPoint.x = 10.9188262;  // x
    iPoint.y = 32.0480347;  // y
    iPoint.z = 22.2950001;  // z
    targetCloud.pts.push_back(iPoint);
    // Element {3}
    iPoint.atomID = 3;      // iatom
    iPoint.x = 7.7338257;   // x
    iPoint.y = 37.564621;   // y
    iPoint.z = 22.2950001;  // z
    targetCloud.pts.push_back(iPoint);
    // Element {4}
    iPoint.atomID = 4;      // iatom
    iPoint.x = 8.1605368;   // x
    iPoint.y = 30.4555359;  // y
    iPoint.z = 25.4799996;  // z
    targetCloud.pts.push_back(iPoint);
    // Element {5}
    iPoint.atomID = 5;      // iatom
    iPoint.x = 6.5680371;   // x
    iPoint.y = 33.2138252;  // y
    iPoint.z = 22.2950001;  // z
    targetCloud.pts.push_back(iPoint);
    // Element {6}
    iPoint.atomID = 6;      // iatom
    iPoint.x = 12.0846195;  // x
    iPoint.y = 36.3988266;  // y
    iPoint.z = 22.2950001;  // z
    targetCloud.pts.push_back(iPoint);
    // Element {7}
    iPoint.atomID = 7;      // iatom
    iPoint.x = 10.3361149;  // x
    iPoint.y = 29.8733234;  // y
    iPoint.z = 23.8880006;  // z
    targetCloud.pts.push_back(iPoint);
    // Element {8}
    iPoint.atomID = 8;      // iatom
    iPoint.x = 7.1511145;   // x
    iPoint.y = 35.3899078;  // y
    iPoint.z = 23.8880006;  // z
    targetCloud.pts.push_back(iPoint);
    // Element {9}
    iPoint.atomID = 9;      // iatom
    iPoint.x = 9.9094058;   // x
    iPoint.y = 36.9824066;  // y
    iPoint.z = 20.7029991;  // z
    targetCloud.pts.push_back(iPoint);
    // Element {10}
    iPoint.atomID = 10;     // iatom
    iPoint.x = 5.985323;    // x
    iPoint.y = 31.039114;   // y
    iPoint.z = 23.8880006;  // z
    targetCloud.pts.push_back(iPoint);
    // Element {11}
    iPoint.atomID = 11;     // iatom
    iPoint.x = 11.5019055;  // x
    iPoint.y = 34.2241135;  // y
    iPoint.z = 23.8880006;  // z
    targetCloud.pts.push_back(iPoint);
    // Element {12}
    iPoint.atomID = 12;     // iatom
    iPoint.x = 8.7436142;   // x
    iPoint.y = 32.6316147;  // y
    iPoint.z = 20.7029991;  // z
    targetCloud.pts.push_back(iPoint);
    // Element {13}
    iPoint.atomID = 13;     // iatom
    iPoint.x = 8.7436142;   // x
    iPoint.y = 32.6316147;  // y
    iPoint.z = 27.073;      // z
    targetCloud.pts.push_back(iPoint);
    // Update nop
    targetCloud.nop = targetCloud.pts.size();
    // box lengths
    targetCloud.box.push_back(50);  // x box length
    targetCloud.box.push_back(50);  // y box length
    targetCloud.box.push_back(50);  // z box length
    // Update the unordered map
    for (int iatom = 0; iatom < targetCloud.nop; iatom++) {
      targetCloud.idIndexMap[iatom] = iatom;
    }  // end of filling the map
    //
    //
    // --------------------------
    // GETTING THE TARGET POINT SET
    // Calculate a neighbour list
    nList = nneigh::neighListO(3.5, &targetCloud, 1);
    // Neighbour list by index
    nList = nneigh::neighbourListByIndex(&targetCloud, nList);
    // Find the vector of vector of rings
    rings = primitive::ringNetwork(nList, 6);
    // init the ringType vector
    ringType.resize(rings.size());
    listDDC = ring::findDDC(rings, &ringType, listHC, &cageList);
    // Save the order of the DDC in a vector
    ddcOrder = pntToPnt::relOrderDDC(0, rings, cageList);
    //
    // --------------------------
    // Now get the absolute orientation of the left (candidate/target) system
    // with respect to the right (template/reference) system test
    //
    std::vector<double> quaternionRot;         // quaternion rotation
    double rmsd1, rmsd2;                       // least RMSD
    std::vector<double> rmsdList1, rmsdList2;  // List of RMSD per atom
    double scale;                              // Scale factor
    //
    //
    // Variables for looping through possible permutations
    //
    std::vector<double> currentQuat;      // quaternion rotation
    double currentRmsd;                   // least RMSD
    std::vector<double> currentRmsdList;  // List of RMSD per atom
    double currentScale;
    // absolute orientation using Horn's algorithm between the target and test
    // set
    int index;
    for (int i = 0; i < 6; i++) {
      // Change the order of the target points somehow!
      //
      targetPointSet = pntToPnt::changeDiaCageOrder(&targetCloud, ddcOrder, i);
      // Shape-matching
      absor::hornAbsOrientation(refPnts, targetPointSet, &currentQuat,
                                &currentRmsd, &currentRmsdList, &currentScale);
      if (i == 0) {
        quaternionRot = currentQuat;
        rmsd1 = currentRmsd;
        rmsdList1 = currentRmsdList;
        scale = currentScale;
        index = 0;
      } else {
        if (currentRmsd < rmsd1) {
          quaternionRot = currentQuat;
          rmsd1 = currentRmsd;
          rmsdList1 = currentRmsdList;
          scale = currentScale;
          index = i;
        }  // update
      }    // Update if this is a better match
    }      // Loop through possible permutations
    // ---------
    std::vector<double>
        selfQuatRot;   // quaternion for the reference set and itself
    double selfScale;  // Scale for the reference set and itself

    // Shape-matching
    absor::hornAbsOrientation(refPnts, refPnts, &selfQuatRot, &rmsd2,
                              &rmsdList2, &selfScale);

    //
    double angDist = gen::angDistDegQuaternions(selfQuatRot, quaternionRot);
    //
    REQUIRE_THAT(angDist, Catch::Matchers::Floating::WithinAbsMatcher(
                              30.0, 0.01));  // Evaluate condition
    // --------------------------
    // --------------------------
  }  // End of given
} // End of scenario
