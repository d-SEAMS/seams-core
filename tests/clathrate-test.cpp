// Internal
#include <absOrientation.hpp>
#include <bulkTUM.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <pntCorrespondence.hpp>
#include <ring.hpp>
#include <topo_bulk.hpp>
#include <bulkClathrate.hpp>
#include <numeric>

// Standard
#include <iostream>

#include <catch2/catch.hpp>
#include <rang.hpp>

SCENARIO("Test the shape-matching of a perfect 5^12 6^4 clathrate cage rotated by 30 degrees",
         "[match]") {
  GIVEN("A rotated clathrate cage (candidate structure)") {
    // --------------------------
    // GETTING THE REFERENCE/TEMPLATE POINT SET
    int dim = 3;   // Number of dimensions
    int nOxy = 28; // Number of O atoms in the clathrate hexadecahedron
    int nop = 84;  // Number of O and H atoms in the clathrate hexadecahedron
    int oxygenAtomType = 1; // LAMMPS type ID for O atoms in the reference structure 
    //
    // Row major (reference point sets)
    Eigen::MatrixXdRowMajor refPntsO(nOxy, dim); // Reference point set of just O atoms (Eigen matrix)
    Eigen::MatrixXdRowMajor refPntsWat(nop, dim); // Reference point set of O H H water atoms (Eigen matrix)
    std::string filePathDumpOonly = "../templates/s2-clathrate/single-cage-O-only.lammpstrj";
    std::string filePathDumpWat = "../templates/s2-clathrate/single-cage-watOnly.lammpstrj";
    // -------------------------------------
    //
    // REFERENCE POINT SET
    // This is row-ordered 
    std::tie(refPntsO, refPntsWat) = clath::buildRefS2CageLammpsTrj(filePathDumpWat, filePathDumpOonly, oxygenAtomType);
    //
    // -------------------------------------
    // GETTING THE ROTATED POINT SET
    molSys::PointCloud<molSys::Point<double>, double>
        targetCloud;  // pointCloud
    molSys::Point<double> iPoint;
    // Column ordered Eigen matrix for the rotated O atom point set 
    Eigen::MatrixXd targetPointSet(nOxy, dim);
    // Fill the pointCloud
    //
    iPoint.type = 1;  // Same for all the points here
    // Element {0}
    iPoint.atomID = 41;         // iatom
    iPoint.molID = 9; // molID 
    iPoint.x = 7.5937367291;  // x
    iPoint.y = 4.9288586199;   // y
    iPoint.z = 6.0233160444;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {1}
    iPoint.atomID = 44;         // iatom
    iPoint.molID = 10; // molID
    iPoint.x = 7.3633886336;  // x
    iPoint.y = 4.1162401801;   // y
    iPoint.z = 8.6870455293;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {2}
    iPoint.atomID = 47;         // iatom
    iPoint.molID = 11; // molID
    iPoint.x = 9.996572126;  // x
    iPoint.y = 4.3055809253;   // y
    iPoint.z = 9.6457430081;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {3}
    iPoint.atomID = 50;         // iatom
    iPoint.molID = 12; // molID
    iPoint.x = 11.7614432412;  // x
    iPoint.y = 5.2981158308;   // y
    iPoint.z = 7.6933335238;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {4}
    iPoint.atomID = 53;         // iatom
    iPoint.molID = 13; // molID
    iPoint.x = 10.2044652858;  // x
    iPoint.y = 5.6949980917;   // y
    iPoint.z = 5.3825971533;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {5}
    iPoint.atomID = 56;         // iatom
    iPoint.molID = 14; // molID
    iPoint.x = 5.2421552804;  // x
    iPoint.y = 5.6518334658;   // y
    iPoint.z = 9.6984569456;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {6}
    iPoint.atomID = 59;         // iatom
    iPoint.molID = 15; // molID
    iPoint.x = 4.2542007054;  // x
    iPoint.y = 7.2846201393;   // y
    iPoint.z = 7.6175500346;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {7}
    iPoint.atomID = 62;         // iatom
    iPoint.molID = 16; // molID
    iPoint.x = 5.8200109197;  // x
    iPoint.y = 6.9701057429;   // y
    iPoint.z = 5.2988053265;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {8}
    iPoint.atomID = 65;         // iatom
    iPoint.molID = 17; // molID
    iPoint.x = 7.3537715067;  // x
    iPoint.y = 9.0759236821;   // y
    iPoint.z = 4.2516259922;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {9}
    iPoint.atomID = 68;         // iatom
    iPoint.molID = 18; // molID
    iPoint.x = 10.0343913562;  // x
    iPoint.y = 8.235249558;   // y
    iPoint.z = 4.2091686876;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {10}
    iPoint.atomID = 71;        // iatom
    iPoint.molID = 19; // molID
    iPoint.x = 9.6785876561;  // x
    iPoint.y = 12.3485536092;   // y
    iPoint.z = 5.8174538112;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {11}
    iPoint.atomID = 74;        // iatom
    iPoint.molID = 20; // molID
    iPoint.x = 9.9499754667;  // x
    iPoint.y = 13.1695490027;   // y
    iPoint.z = 8.4516453084;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {12}
    iPoint.atomID = 77;        // iatom
    iPoint.molID = 21; // molID
    iPoint.x = 7.4033636906;  // x
    iPoint.y = 12.8878722476;   // y
    iPoint.z = 9.6442399292;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {13}
    iPoint.atomID = 80;        // iatom
    iPoint.molID = 22; // molID
    iPoint.x = 5.5594338226;  // x
    iPoint.y = 11.9941133355;   // y
    iPoint.z = 7.6892320116;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {14}
    iPoint.atomID = 83;        // iatom
    iPoint.molID = 23; // molID
    iPoint.x = 7.0304681436;  // x
    iPoint.y = 11.6579787992;   // y
    iPoint.z = 5.3331503037;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {15}
    iPoint.atomID = 86;        // iatom
    iPoint.molID = 24; // molID
    iPoint.x = 12.0636261256;  // x
    iPoint.y = 11.6986624899 ;   // y
    iPoint.z = 9.6001962219;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {16}
    iPoint.atomID = 89;        // iatom
    iPoint.molID = 25; // molID
    iPoint.x = 12.9869681953;  // x
    iPoint.y = 9.9725925792;   // y
    iPoint.z = 7.5755747819;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {17}
    iPoint.atomID = 92;        // iatom
    iPoint.molID = 26; // molID
    iPoint.x = 11.5343919355;  // x
    iPoint.y = 10.3673111086;   // y
    iPoint.z = 5.2299192572;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {18}
    iPoint.atomID = 95;        // iatom
    iPoint.molID = 27; // molID
    iPoint.x = 5.0098242829;  // x
    iPoint.y = 9.6456766276;   // y
    iPoint.z = 11.3443477392;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {19}
    iPoint.atomID = 98;        // iatom
    iPoint.molID = 28; // molID
    iPoint.x = 4.0682759589;  // x
    iPoint.y = 9.8559840063;   // y
    iPoint.z = 8.720098333;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {20}
    iPoint.atomID = 101;        // iatom
    iPoint.molID = 29; // molID
    iPoint.x = 6.9914005664;  // x
    iPoint.y = 11.4847835543;   // y
    iPoint.z = 12.0273021347;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {21}
    iPoint.atomID = 104;        // iatom
    iPoint.molID = 30; // molID
    iPoint.x = 5.8113593439;  // x
    iPoint.y = 7.0652181081;   // y
    iPoint.z = 12.0425269542;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {22}
    iPoint.atomID = 107;        // iatom
    iPoint.molID = 31; // molID
    iPoint.x = 8.3735199968;  // x
    iPoint.y = 7.3990394087;   // y
    iPoint.z = 13.151513935;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {23}
    iPoint.atomID = 110;        // iatom
    iPoint.molID = 32; // molID
    iPoint.x = 9.1561932071;  // x
    iPoint.y = 10.0777760966;   // y
    iPoint.z = 13.1118134976;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {24}
    iPoint.atomID = 113;        // iatom
    iPoint.molID = 33; // molID
    iPoint.x = 12.4670069217;  // x
    iPoint.y = 7.570750914;   // y
    iPoint.z = 11.3668569684;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {25}
    iPoint.atomID = 116;        // iatom
    iPoint.molID = 34; // molID
    iPoint.x = 13.2967536635;  // x
    iPoint.y = 7.4183025368;   // y
    iPoint.z = 8.7082803336;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {26}
    iPoint.atomID = 119;        // iatom
    iPoint.molID = 35; // molID
    iPoint.x = 10.3525245971;  // x
    iPoint.y = 5.8336554191;   // y
    iPoint.z = 11.9513600459;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {27}
    iPoint.atomID = 122;        // iatom
    iPoint.molID = 36; // molID
    iPoint.x = 11.7014756194;  // x
    iPoint.y = 10.175494955;   // y
    iPoint.z = 11.9611383524;     // z
    targetCloud.pts.push_back(iPoint);
    // Update nop
    targetCloud.nop = targetCloud.pts.size();
    // box lengths
    targetCloud.box.push_back(17.364);  // x box length
    targetCloud.box.push_back(17.364);  // y box length
    targetCloud.box.push_back(17.364);  // z box length
    // Update the unordered map
    for (int iatom = 0; iatom < targetCloud.nop; iatom++) {
      targetCloud.idIndexMap[targetCloud.pts[iatom].atomID] = iatom;
    }  // end of filling the map
    //
    //
    // --------------------------
    // GETTING THE TARGET POINT SET
    // Column major 
    Eigen::MatrixXdRowMajor targetPnts(nOxy, dim); // Target point set of just O atoms (Eigen matrix)
    std::vector<int> atomIndexList; // The vector with the atom indices 
                                    // corresponding to the order of the atoms in the Eigen target point set 
    
    // Get the initial ordering of atomIndexList
    // atomIndexList should be in ascending order
    for (int i = 0; i < nOxy; ++i)
    {
        atomIndexList.push_back(i); // index from 0 to nOxy-1
    } // this is automatically in ascending order 

    // We need to loop through all the permutations of atomIndexList
    // and get the absolute orientation of the left (candidate/target) system
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

    // Get all the possible permutations, starting with the sorted
    // atomIndexList vector corresponding to atom indices in the targetCloud
    int count=0; // for looping through all permutations
    int permutNum; // permutation number, starting from 0  

    do {
        // reordered atomIndexList
        // Create the target Eigen matrix point set 
        targetPnts = pntToPnt::fillTargetEigenPointSet(targetCloud,
            atomIndexList, nOxy, dim); 
        //
        // Shape-matching
        absor::hornAbsOrientationRowMajor(refPntsO, targetPnts, &currentQuat,
            &currentRmsd, &currentRmsdList, &currentScale); 

        // Update if the new permutation is a better match 
        if (count == 0) {
        quaternionRot = currentQuat;
        rmsd1 = currentRmsd;
        rmsdList1 = currentRmsdList;
        scale = currentScale;
        permutNum = count; 
        break;
        } else {
            if (currentRmsd < rmsd1) {
          quaternionRot = currentQuat;
          rmsd1 = currentRmsd;
          rmsdList1 = currentRmsdList;
          scale = currentScale;
          permutNum = count; 
            }  // update
        }    // Update if this is a better match

        // Update count 
        count++; 
    } while (std::next_permutation(atomIndexList.begin(), atomIndexList.end()));

    // ---------
    std::vector<double>
        selfQuatRot;   // quaternion for the reference set and itself
    double selfScale;  // Scale for the reference set and itself

    // Shape-matching
    absor::hornAbsOrientationRowMajor(refPntsO, refPntsO, &selfQuatRot, &rmsd2,
                              &rmsdList2, &selfScale);

    //
    double angDist = gen::angDistDegQuaternions(selfQuatRot, quaternionRot);
    //
    REQUIRE_THAT(angDist, Catch::Matchers::Floating::WithinAbsMatcher(
                              30.0, 0.01));  // Evaluate condition
    // --------------------------
  }  // End of given
}  // End of scenario

SCENARIO("Test the shape-matching of a perfect 5^12 6^4 clathrate cage with one rotated by 30 degrees, using primitive rings",
         "[match]") {
  GIVEN("A rotated clathrate cage (candidate structure)") {
    // --------------------------
    // GETTING THE REFERENCE/TEMPLATE POINT SET
    int dim = 3;   // Number of dimensions
    int nOxy = 28; // Number of O atoms in the clathrate hexadecahedron
    int oxygenAtomType = 1; // LAMMPS type ID for O atoms in the reference structure 
    //
    // Ring stuff 
    std::vector<std::vector<int>> nList;  // Neighbour list
    std::vector<std::vector<int>> rings, ringsRef;  // Rings
    std::vector<std::vector<int>>
      ringsHex, ringsAll;    // Vector of vectors of rings
    //
    // Row major (reference point sets)
    Eigen::MatrixXdRowMajor refPntsO(nOxy, dim); // Reference point set of just O atoms (Eigen matrix)
    // Row ordered Eigen matrix for the rotated O atom point set 
    Eigen::MatrixXd targetPointSet(nOxy, dim);
    // PointClouds 
    molSys::PointCloud<molSys::Point<double>, double>
        refCloud;  // pointCloud
    molSys::PointCloud<molSys::Point<double>, double>
        targetCloud;  // pointCloud for the target point set 
    std::string filePathDumpOonly = "../templates/s2-clathrate/single-cage-O-only.lammpstrj";
    // RMSD stuff 
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
    // -------------------------------------
    //
    // REFERENCE POINT SET
    // Get the reference PointCloud 
    refCloud = sinp::readLammpsTrjO(filePathDumpOonly, 1, &refCloud, oxygenAtomType);
    // This is row-ordered 
    std::tie(refCloud, ringsRef, refPntsO) = clath::buildRefS2Cage(filePathDumpOonly, oxygenAtomType);
    //
    // --------------------------
    // GETTING THE TARGET POINT SET
    molSys::Point<double> iPoint;
    // Fill the pointCloud
    //
    iPoint.type = 1;  // Same for all the points here
    // Element {0}
    iPoint.atomID = 41;         // iatom
    iPoint.molID = 9; // molID 
    iPoint.x = 7.5937367291;  // x
    iPoint.y = 4.9288586199;   // y
    iPoint.z = 6.0233160444;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {1}
    iPoint.atomID = 44;         // iatom
    iPoint.molID = 10; // molID
    iPoint.x = 7.3633886336;  // x
    iPoint.y = 4.1162401801;   // y
    iPoint.z = 8.6870455293;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {2}
    iPoint.atomID = 47;         // iatom
    iPoint.molID = 11; // molID
    iPoint.x = 9.996572126;  // x
    iPoint.y = 4.3055809253;   // y
    iPoint.z = 9.6457430081;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {3}
    iPoint.atomID = 50;         // iatom
    iPoint.molID = 12; // molID
    iPoint.x = 11.7614432412;  // x
    iPoint.y = 5.2981158308;   // y
    iPoint.z = 7.6933335238;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {4}
    iPoint.atomID = 53;         // iatom
    iPoint.molID = 13; // molID
    iPoint.x = 10.2044652858;  // x
    iPoint.y = 5.6949980917;   // y
    iPoint.z = 5.3825971533;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {5}
    iPoint.atomID = 56;         // iatom
    iPoint.molID = 14; // molID
    iPoint.x = 5.2421552804;  // x
    iPoint.y = 5.6518334658;   // y
    iPoint.z = 9.6984569456;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {6}
    iPoint.atomID = 59;         // iatom
    iPoint.molID = 15; // molID
    iPoint.x = 4.2542007054;  // x
    iPoint.y = 7.2846201393;   // y
    iPoint.z = 7.6175500346;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {7}
    iPoint.atomID = 62;         // iatom
    iPoint.molID = 16; // molID
    iPoint.x = 5.8200109197;  // x
    iPoint.y = 6.9701057429;   // y
    iPoint.z = 5.2988053265;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {8}
    iPoint.atomID = 65;         // iatom
    iPoint.molID = 17; // molID
    iPoint.x = 7.3537715067;  // x
    iPoint.y = 9.0759236821;   // y
    iPoint.z = 4.2516259922;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {9}
    iPoint.atomID = 68;         // iatom
    iPoint.molID = 18; // molID
    iPoint.x = 10.0343913562;  // x
    iPoint.y = 8.235249558;   // y
    iPoint.z = 4.2091686876;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {10}
    iPoint.atomID = 71;        // iatom
    iPoint.molID = 19; // molID
    iPoint.x = 9.6785876561;  // x
    iPoint.y = 12.3485536092;   // y
    iPoint.z = 5.8174538112;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {11}
    iPoint.atomID = 74;        // iatom
    iPoint.molID = 20; // molID
    iPoint.x = 9.9499754667;  // x
    iPoint.y = 13.1695490027;   // y
    iPoint.z = 8.4516453084;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {12}
    iPoint.atomID = 77;        // iatom
    iPoint.molID = 21; // molID
    iPoint.x = 7.4033636906;  // x
    iPoint.y = 12.8878722476;   // y
    iPoint.z = 9.6442399292;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {13}
    iPoint.atomID = 80;        // iatom
    iPoint.molID = 22; // molID
    iPoint.x = 5.5594338226;  // x
    iPoint.y = 11.9941133355;   // y
    iPoint.z = 7.6892320116;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {14}
    iPoint.atomID = 83;        // iatom
    iPoint.molID = 23; // molID
    iPoint.x = 7.0304681436;  // x
    iPoint.y = 11.6579787992;   // y
    iPoint.z = 5.3331503037;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {15}
    iPoint.atomID = 86;        // iatom
    iPoint.molID = 24; // molID
    iPoint.x = 12.0636261256;  // x
    iPoint.y = 11.6986624899 ;   // y
    iPoint.z = 9.6001962219;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {16}
    iPoint.atomID = 89;        // iatom
    iPoint.molID = 25; // molID
    iPoint.x = 12.9869681953;  // x
    iPoint.y = 9.9725925792;   // y
    iPoint.z = 7.5755747819;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {17}
    iPoint.atomID = 92;        // iatom
    iPoint.molID = 26; // molID
    iPoint.x = 11.5343919355;  // x
    iPoint.y = 10.3673111086;   // y
    iPoint.z = 5.2299192572;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {18}
    iPoint.atomID = 95;        // iatom
    iPoint.molID = 27; // molID
    iPoint.x = 5.0098242829;  // x
    iPoint.y = 9.6456766276;   // y
    iPoint.z = 11.3443477392;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {19}
    iPoint.atomID = 98;        // iatom
    iPoint.molID = 28; // molID
    iPoint.x = 4.0682759589;  // x
    iPoint.y = 9.8559840063;   // y
    iPoint.z = 8.720098333;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {20}
    iPoint.atomID = 101;        // iatom
    iPoint.molID = 29; // molID
    iPoint.x = 6.9914005664;  // x
    iPoint.y = 11.4847835543;   // y
    iPoint.z = 12.0273021347;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {21}
    iPoint.atomID = 104;        // iatom
    iPoint.molID = 30; // molID
    iPoint.x = 5.8113593439;  // x
    iPoint.y = 7.0652181081;   // y
    iPoint.z = 12.0425269542;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {22}
    iPoint.atomID = 107;        // iatom
    iPoint.molID = 31; // molID
    iPoint.x = 8.3735199968;  // x
    iPoint.y = 7.3990394087;   // y
    iPoint.z = 13.151513935;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {23}
    iPoint.atomID = 110;        // iatom
    iPoint.molID = 32; // molID
    iPoint.x = 9.1561932071;  // x
    iPoint.y = 10.0777760966;   // y
    iPoint.z = 13.1118134976;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {24}
    iPoint.atomID = 113;        // iatom
    iPoint.molID = 33; // molID
    iPoint.x = 12.4670069217;  // x
    iPoint.y = 7.570750914;   // y
    iPoint.z = 11.3668569684;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {25}
    iPoint.atomID = 116;        // iatom
    iPoint.molID = 34; // molID
    iPoint.x = 13.2967536635;  // x
    iPoint.y = 7.4183025368;   // y
    iPoint.z = 8.7082803336;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {26}
    iPoint.atomID = 119;        // iatom
    iPoint.molID = 35; // molID
    iPoint.x = 10.3525245971;  // x
    iPoint.y = 5.8336554191;   // y
    iPoint.z = 11.9513600459;     // z
    targetCloud.pts.push_back(iPoint);
    // Element {27}
    iPoint.atomID = 122;        // iatom
    iPoint.molID = 36; // molID
    iPoint.x = 11.7014756194;  // x
    iPoint.y = 10.175494955;   // y
    iPoint.z = 11.9611383524;     // z
    targetCloud.pts.push_back(iPoint);
    // Update nop
    targetCloud.nop = targetCloud.pts.size();
    // box lengths
    targetCloud.box.push_back(17.364);  // x box length
    targetCloud.box.push_back(17.364);  // y box length
    targetCloud.box.push_back(17.364);  // z box length
    // Update the unordered map
    for (int iatom = 0; iatom < targetCloud.nop; iatom++) {
      targetCloud.idIndexMap[targetCloud.pts[iatom].atomID] = iatom;
    }  // end of filling the map
    //
    // -----------------------
    // TARGET POINT SET 
    // For the target point set 
    // Calculate a neighbour list
    nList = nneigh::neighListO(3.5, &targetCloud, 1);
    // Neighbour list by index
    nList = nneigh::neighbourListByIndex(&targetCloud, nList);
    // Find the vector of vector of rings
    ringsAll = primitive::ringNetwork(nList, 6); 
    // Get just the 6-membered rings 
    ringsHex = ring::getSingleRingSize(ringsAll, 6);
    // SHAPE-MATCHING 
    // 6-membered rings 
    std::vector<int> tempRing, revRing;
    std::vector<int> currentLastRing;

    // Loop through the rings 
    for (auto& iring : ringsHex)
    {
        tempRing = iring;
        revRing = iring;
        // Reverse the ring
        std::reverse(revRing.begin(), revRing.end());
        // No reversal 
        for (int i = 0; i < ringsHex[0].size(); ++i)
        {
            if (i==0)
            {
                // -------------------
                clath::matchClathrateLastRing(rings, tempRing, ringsRef, 
                    targetCloud, refCloud,&currentQuat,
                    &currentRmsd, &currentRmsdList, &currentScale);
                // Update for the first time 
                quaternionRot = currentQuat;
                rmsd1 = currentRmsd;
                rmsdList1 = currentRmsdList;
                scale = currentScale;
                currentLastRing = tempRing;
                // -------------------
            } // first step 
            else{
                // Change the order of the ring 
                rotate(tempRing.begin(), tempRing.begin()+1, tempRing.end());
                // Shape-matching 
                clath::matchClathrateLastRing(rings, tempRing, ringsRef, 
                    targetCloud, refCloud,&currentQuat,
                    &currentRmsd, &currentRmsdList, &currentScale);
                // Update if currentRmsd is less than rmsd1
                if (currentRmsd < rmsd1)
                {
                    quaternionRot = currentQuat;
                    rmsd1 = currentRmsd;
                    rmsdList1 = currentRmsdList;
                    scale = currentScale;
                    currentLastRing = tempRing;
                } // end of update 
            } // all other steps apart from the first 
        } // go through 12 arrangements of the ring 

        // Reversed ring  
        for (int i = 0; i < ringsHex[0].size(); ++i)
        {
            if (i==0)
            {
                // -------------------
                clath::matchClathrateLastRing(rings, revRing, ringsRef, 
                    targetCloud, refCloud,&currentQuat,
                    &currentRmsd, &currentRmsdList, &currentScale);
                if (currentRmsd < rmsd1)
                {
                    quaternionRot = currentQuat;
                    rmsd1 = currentRmsd;
                    rmsdList1 = currentRmsdList;
                    scale = currentScale;
                    currentLastRing = revRing;
                } // end of update 
                // -------------------
            } // first step 
            else{
                // Change the order of the ring 
                rotate(revRing.begin(), revRing.begin()+1, revRing.end());
                // Shape-matching 
                clath::matchClathrateLastRing(rings, revRing, ringsRef, 
                    targetCloud, refCloud,&currentQuat,
                    &currentRmsd, &currentRmsdList, &currentScale);
                // Update if currentRmsd is less than rmsd1
                if (currentRmsd < rmsd1)
                {
                    quaternionRot = currentQuat;
                    rmsd1 = currentRmsd;
                    rmsdList1 = currentRmsdList;
                    scale = currentScale;
                    currentLastRing = revRing;
                } // end of update 
            } // all other steps apart from the first 
        } // go through 12 arrangements of the ring 

        rings.push_back(currentLastRing);

    } // end of loop through ringsHex
    // ---------
    // Get the last 4 elements, corresponding to the 5th ring in ringsRef
    // ringsRef[4]

    // 28 atom indices in targetCloud
    std::vector<int> targetAtomIndices; 

    // Get all the 28 indices of the atoms in targetCloud
    // This could be more complicated in a large system, where
    // the closest 28 water molecules to the THF center of mass might 
    // be chosen 
    for (int iatom = 0; iatom < targetCloud.nop; iatom++)
    {
        targetAtomIndices.push_back(iatom);
    } // overkill!

    // Flattened rings vector, which will contain all the indices
    // in the rings vector for the targetCloud
    auto flattenedRingsVec = std::accumulate(rings.begin(), rings.end(), decltype(rings)::value_type{},
            [](auto &x, auto &y) {
        x.insert(x.end(), y.begin(), y.end());
        return x;
    });

    // Sort the flattened rings vector and targetAtomIndices
    std::sort(targetAtomIndices.begin(), targetAtomIndices.end());
    std::sort(flattenedRingsVec.begin(), flattenedRingsVec.end());

    std::vector<int> lastTargetRing; // The last ring with 4 elements

    // Elements not in common 
    std::set_symmetric_difference(
        targetAtomIndices.begin(), targetAtomIndices.end(),
        flattenedRingsVec.begin(), flattenedRingsVec.end(),
        std::back_inserter(lastTargetRing));

    // Get all possible permutations and do shape-matching
    // of the last ring 
    std::sort(lastTargetRing); // should be in ascending order
    int count=0; // for looping through all permutations

    // Go through all the permutations 
    do {
        // reordered lastTargetRing
        if (count==0)
        {
            // -------------------
            clath::matchClathrateLastRing(rings, lastTargetRing, ringsRef, 
                targetCloud, refCloud,&currentQuat,
                &currentRmsd, &currentRmsdList, &currentScale);
            // Update for the first time 
            quaternionRot = currentQuat;
            rmsd1 = currentRmsd;
            rmsdList1 = currentRmsdList;
            scale = currentScale;
            currentLastRing = lastTargetRing;
            // -------------------
        } // first permutation
        else{
            // Shape-matching 
            clath::matchClathrateLastRing(rings, lastTargetRing, ringsRef, 
                targetCloud, refCloud,&currentQuat,
                &currentRmsd, &currentRmsdList, &currentScale);
            // Update if currentRmsd is less than rmsd1
            if (currentRmsd < rmsd1)
            {
                quaternionRot = currentQuat;
                rmsd1 = currentRmsd;
                rmsdList1 = currentRmsdList;
                scale = currentScale;
                currentLastRing = lastTargetRing;
            } // end of update
        } // subsequent permutations after the first  

        count++; // update the number of times the loop has run 
    } while (std::next_permutation(lastTargetRing.begin(), lastTargetRing.end()));

    // Presumably currentLastRing is the best match 
    rings.push_back(currentLastRing);

    // BREAK HERE AND TEST 
    // -------------------------
    std::vector<double>
        selfQuatRot;   // quaternion for the reference set and itself
    double selfScale;  // Scale for the reference set and itself

    // Shape-matching
    absor::hornAbsOrientationRowMajor(refPntsO, refPntsO, &selfQuatRot, &rmsd2,
                              &rmsdList2, &selfScale);

    //
    double angDist = gen::angDistDegQuaternions(selfQuatRot, selfQuatRot);
    //
    REQUIRE_THAT(angDist, Catch::Matchers::Floating::WithinAbsMatcher(
                              0.0, 0.01));  // Evaluate condition
    // --------------------------
  }  // End of given
}  // End of scenario