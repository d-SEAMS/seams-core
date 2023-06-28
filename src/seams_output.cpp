//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#include <seams_input.hpp>
#include <seams_output.hpp>

/**
 * @details  Prints out a LAMMPS data file, with some default options. Only
 * Oxygen atoms are printed out
 */
int sout::writeLAMMPSdata(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> rings, std::vector<std::vector<int>> bonds,
    std::string filename) {
  std::ofstream outputFile;
  std::vector<int> atoms;         // Holds all atom IDs to print
  int ringSize = rings[0].size(); // Ring size of each ring in rings
  int iatom;                      // Index, not atom ID
  bool padAtoms = false;          // Add extra atoms if the atom IDs are skipped
  int prevAtomID = 0;             // Check for previous atom ID
  int dummyAtoms = 0;             // Number of dummy atoms to fill
  int dummyID;
  int jatom; // Array index is 1 less than the ID (index for dummy atom)
  // ----------------
  // Return if there are no rings
  if (rings.size() == 0) {
    return 1;
  }
  // ----------------
  // Otherwise create file
  // Create output dir if it doesn't exist already
  const char *path = "../output"; // relative to the build directory
  fs::path dir(path);
  // if (fs::create_directory(dir)) {
  //   std::cerr << "Output directory created\n";
  // }
  // ----------------
  // Get all the unique atomIDs of the atoms in the rings of this type
  // Put all atom IDs into one 1-D vector
  size_t total_size{0};
  // Get the total number of atoms (repeated)
  total_size = rings.size() * ringSize;
  // Reserve this size inside atoms
  atoms.reserve(total_size);
  // Fill up all these atom IDs
  for (int iring = 0; iring < rings.size(); iring++) {
    std::move(rings[iring].begin(), rings[iring].end(),
              std::back_inserter(atoms));
  } // end of loop through all rings in the list

  // Sort the array according to atom ID, which will be needed to get the
  // unique IDs and to remove duplicates
  sort(atoms.begin(), atoms.end());
  // Get the unique atom IDs
  auto ip = std::unique(atoms.begin(), atoms.end());
  // Resize the vector to remove undefined terms
  atoms.resize(std::distance(atoms.begin(), ip));
  // If the number of atoms is less than the total nop, add dummy atoms
  if (atoms.size() != yCloud->nop) {
    padAtoms = true;
  }
  // ----------------
  // Write output to file inside the output directory
  outputFile.open("../output/" + filename);
  // FORMAT:
  //  Comment Line
  //  4 atoms
  //  4 bonds
  //  0 angles
  //  0 dihedrals
  //  0 impropers
  //  1 atom types
  //  1 bond types
  //  0 angle types
  //  0 dihedral types
  //  0 improper types
  //  -1.124000 52.845002  xlo xhi
  //  0.000000 54.528999  ylo yhi
  //  1.830501 53.087501  zlo zhi

  //  Masses

  //  1 15.999400 # O

  //  Atoms

  // 1 1 1 0 20.239  1.298 6.873 # O
  // 2 1 1 0 0 5.193 6.873 # O
  // 3 1 1 0 2.249 1.298 6.873 # O

  // -------
  // Write the header
  // Write comment line
  outputFile << "Written out by D-SEAMS\n";
  // Write out the number of atoms
  outputFile << atoms[atoms.size() - 1] << " "
             << "atoms"
             << "\n";
  // Number of bonds
  outputFile << bonds.size() << " bonds"
             << "\n";
  outputFile << "0 angles\n0 dihedrals\n0 impropers\n";
  // If padded atoms are required, two atom types will be required
  if (padAtoms) {
    outputFile << "2 atom types\n";
  } else {
    outputFile << "1 atom types\n";
  } // end of atom types
  outputFile
      << "1 bond types\n0 angle types\n0 dihedral types\n0 improper types\n";
  // Box lengths
  outputFile << "0 " << yCloud->box[0] << " xlo xhi\n";
  outputFile << "0 " << yCloud->box[1] << " ylo yhi\n";
  outputFile << "0 " << yCloud->box[2] << " zlo zhi\n";
  // Masses
  outputFile << "\nMasses\n\n";
  outputFile << "1 15.999400 # O\n";
  if (padAtoms) {
    outputFile << "2 1.0 # dummy\n";
  }
  // Atoms
  outputFile << "\nAtoms\n\n";
  // -------
  // Write out the atom coordinates
  // Loop through atoms
  for (int i = 0; i < atoms.size(); i++) {
    iatom = atoms[i] - 1; // The actual index is one less than the ID
    // -----------
    // Pad out
    // Fill in dummy atoms if some have been skipped
    if (atoms[i] != prevAtomID + 1) {
      dummyAtoms = atoms[i] - prevAtomID - 1;
      dummyID = prevAtomID;
      // Loop to write out dummy atoms
      for (int j = 0; j < dummyAtoms; j++) {
        dummyID++;
        jatom = dummyID - 1;
        // 1 molecule-tag atom-type q x y z
        outputFile << dummyID << " " << yCloud->pts[jatom].molID << " 2 0 "
                   << yCloud->pts[jatom].x << " " << yCloud->pts[jatom].y << " "
                   << yCloud->pts[jatom].z << "\n";
      } // end of dummy atom write-out
    }   // end of check for dummy atom printing
    // -----------
    // Write out coordinates
    // 1 molecule-tag atom-type q x y z
    outputFile << atoms[i] << " " << yCloud->pts[iatom].molID << " 1 0 "
               << yCloud->pts[iatom].x << " " << yCloud->pts[iatom].y << " "
               << yCloud->pts[iatom].z << "\n";
    // update the previous atom ID
    prevAtomID = atoms[i];
  } // end of loop through all atoms in atomID

  // Print the bonds now!
  outputFile << "\nBonds\n\n";
  // Loop through all bonds
  for (int ibond = 0; ibond < bonds.size(); ibond++) {
    outputFile << ibond + 1 << " 1 " << bonds[ibond][0] << " "
               << bonds[ibond][1] << "\n";
  }

  // Once the datafile has been printed, exit
  return 0;
}

/**
 * @details Function for printing out ring info, when there is no
 *  volume slice
 *  Uses Boost!
 */
int sout::writeRings(std::vector<std::vector<int>> rings,
                     std::string filename) {
  std::ofstream outputFile;
  // ----------------
  // Create output dir if it doesn't exist already
  const char *path = "../output"; // relative to the build directory
  fs::path dir(path);
  // if (fs::create_directory(dir)) {
  //   std::cerr << "Output directory created\n";
  // }
  // ----------------
  // Write output to file inside the output directory
  outputFile.open("../output/" + filename);

  // Format:
  // 272    214    906   1361    388      1

  for (int iring = 0; iring < rings.size(); iring++) {
    // Otherwise, write out to the file
    for (int k = 0; k < rings[iring].size(); k++) {
      outputFile << rings[iring][k] << " ";
    } // end of loop through ring elements
    outputFile << "\n";
  } // end of loop through rings

  outputFile.close();

  return 0;
}

/**
 * @details Function for printing out each info, when there is no volume slice
 * Uses Boost!
 */
int sout::writePrisms(
    std::vector<int> *basal1, std::vector<int> *basal2, int prismNum,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud) {
  std::ofstream outputFile;
  std::string number = std::to_string(prismNum);
  std::string filename = "prism" + number + ".dat";
  int ringSize =
      (*basal1).size(); // Size of the ring; each ring contains n elements
  int iatomIndex;       // index of atom coordinate being written out

  // ----------------
  // Create output dir if it doesn't exist already
  const char *path = "../output/prisms"; // relative to the build directory
  fs::path dir(path);
  // if (fs::create_directory(dir)) {
  //   std::cerr << "Output Prism directory created\n";
  // }
  // ----------------
  // Write output to file inside the output directory
  outputFile.open("../output/prisms/" + filename);

  // Format:
  // x y z coordinates of each node

  // For basal 1
  for (int iring = 0; iring < ringSize; iring++) {
    iatomIndex = (*basal1)[iring]; // C++ indices are one less
    // Write the coordinates out to the file
    outputFile << yCloud->pts[iatomIndex].x << " ";
    outputFile << yCloud->pts[iatomIndex].y << " ";
    outputFile << yCloud->pts[iatomIndex].z << " ";

    outputFile << "\n";
  } // end of loop through basal1

  // For basal 2
  for (int iring = 0; iring < ringSize; iring++) {
    iatomIndex = (*basal2)[iring]; // C++ indices are one less
    // Write the coordinates out to the file
    outputFile << yCloud->pts[iatomIndex].x << " ";
    outputFile << yCloud->pts[iatomIndex].y << " ";
    outputFile << yCloud->pts[iatomIndex].z << " ";

    outputFile << "\n";
  } // end of loop through basal1

  outputFile.close();

  // ---- Print out all the coordinates of a single ring, for prismNum=1 only
  if (prismNum == 1) {
    outputFile.open("../output/prisms/singleRing.dat");
    // For basal 1
    for (int iring = 0; iring < ringSize; iring++) {
      iatomIndex = (*basal1)[iring]; // C++ indices are one less
      // Write the coordinates out to the file
      outputFile << yCloud->pts[iatomIndex].x << " ";
      outputFile << yCloud->pts[iatomIndex].y << " ";
      outputFile << yCloud->pts[iatomIndex].z << " ";

      outputFile << "\n";
    } // end of loop through basal1
    outputFile.close();
  }

  return 0;
}

/**
 * @details  Function for writing out all the cages of all types to a folder in
 * the output
 */
int sout::writeAllCages(
    std::string path, std::vector<cage::Cage> *cageList,
    std::vector<std::vector<int>> rings, std::vector<std::vector<int>> nList,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    int currentFrame) {
  int numDDC;                          // Number of DDCs
  int numHC;                           // Number of HCs
  int numMC;                           // Number of MCs
  int totalCages = (*cageList).size(); // Total number of cages
  cage::cageType type;                 // Current cage type

  // ---------------------------------
  // Error handling!
  if (totalCages == 0) {
    std::cerr << "There are no cages to print.\n";
    return 1;
  }
  // ---------------------------------
  // Init
  numDDC = 0;
  numHC = 0;
  numMC = 0;

  // Loop through every cage
  for (int icage = 0; icage < totalCages; icage++) {
    type = (*cageList)[icage].type;
    // ------
    // Add to the cage type and write out to the appropriate folders
    // Hexagonal Cages
    if (type == cage::cageType::HexC) {
      numHC++;
      sout::writeEachCage((*cageList)[icage].rings, numHC, type, rings, yCloud);
      sout::writeBasalRingsHex((*cageList)[icage].rings, numHC, nList, rings);
    } // end of write out of HCs
    // Double diamond Cages
    else if (type == cage::cageType::DoubleDiaC) {
      numDDC++;
      sout::writeEachCage((*cageList)[icage].rings, numDDC, type, rings,
                          yCloud);
    } // end of write out of DDCs
    // // Mixed Cages
    // else if (type == cage::Mixed) {
    //   numMC++;
    //   sout::writeEachCage((*cageList)[icage].rings, numMC, type, rings,
    //   yCloud);
    // }  // end of write out of MCs
    // // Error
    else {
      std::cerr << "The cage is of the wrong type\n";
      continue;
    } // some error
    // ------
  } // end of loop through all cages

  return 0;
}

/**
 * @details Function for printing out each cage's coordinates, when there is no
 * volume slice  Uses Boost!
 */
int sout::writeEachCage(
    std::vector<int> currentCage, int cageNum, cage::cageType type,
    std::vector<std::vector<int>> rings,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud) {
  std::ofstream outputFile;
  std::string number = std::to_string(cageNum);
  std::string filename = "cage" + number + ".dat";
  int ringSize =
      rings[0].size();        // Size of the ring; each ring contains n elements
  int iatomIndex;             // index of atom coordinate being written out
  std::string actualCageType; // is icage a DDC, HC or MC?
  char cageChar[100];         // is icage a DDC, HC or MC?
  int iring;                  // Ring index of the current ring

  if (type == cage::cageType::HexC) {
    strcpy(cageChar, "../output/cages/hexCages");
    actualCageType = "hexCages";
  } else if (type == cage::cageType::DoubleDiaC) {
    strcpy(cageChar, "../output/cages/doubleDiaCages");
    actualCageType = "doubleDiaCages";
  } else {
    // throw error
    std::cerr << "The cage is of the wrong type. Exit\n";
    return 1;
  } //

  // ----------------
  // Create output dir if it doesn't exist already
  const char *path = "../output/cages"; // relative to the build directory
  fs::path dir(path);
  // if (fs::create_directory(dir)) {
  //   std::cerr << "Output Cage directory created\n";
  // }
  // ----------------
  // Subdirectory

  const fs::path path1(cageChar);
  // fs::create_directories(path1);

  // Write output to file inside the output directory
  std::string fileOutNameFull =
      "../output/cages/" + actualCageType + "/" + filename;
  outputFile.open(fileOutNameFull);

  // Format:
  // x y z coordinates of each node

  // For hexagonal cages:
  if (type == cage::cageType::HexC) {
    // Print out only basal ring atoms, since they describe the outer structure
    // The first two rings are basal rings
    for (int i = 0; i < 2; i++) {
      iring = currentCage[i]; // Current iring
      // Get every node of iring
      for (int j = 0; j < ringSize; j++) {
        iatomIndex = rings[iring][j] - 1; // C++ indices are one less
        // Write out the coordinates to the file
        outputFile << yCloud->pts[iatomIndex].x << " ";
        outputFile << yCloud->pts[iatomIndex].y << " ";
        outputFile << yCloud->pts[iatomIndex].z << " ";

        outputFile << "\n";
      } // end of loop through iring
    }   // end of getting every ring in the current cage
  }     // end of printing basal ring atoms for hexagonal cages
  // For currentCage
  else {
    // Loop through all cages (could be duplicates) TODO: remove duplicates
    for (int i = 0; i < currentCage.size(); i++) {
      iring = currentCage[i]; // Current iring
      // Get every node of iring
      for (int j = 0; j < ringSize; j++) {
        iatomIndex = rings[iring][j] - 1; // C++ indices are one less
        // Write out the coordinates to the file
        outputFile << yCloud->pts[iatomIndex].x << " ";
        outputFile << yCloud->pts[iatomIndex].y << " ";
        outputFile << yCloud->pts[iatomIndex].z << " ";

        outputFile << "\n";
      } // end of loop through iring
    }   // end of getting every ring in the current cage
  }     // end of cage printing (has duplicates)

  // Close the output file
  outputFile.close();

  return 0;
}

/**
 * @details Function for printing out the basal rings only of the hexagonal cage
 * described by the number cageNum Uses Boost!
 */
int sout::writeBasalRingsHex(std::vector<int> currentCage, int cageNum,
                             std::vector<std::vector<int>> nList,
                             std::vector<std::vector<int>> rings) {
  std::ofstream outputFile;
  std::string number = std::to_string(cageNum);
  std::string filename = "basalRings" + number + ".dat";
  int ringSize =
      rings[0].size(); // Size of the ring; each ring contains n elements
  int iatomIndex; // index of atom coordinate being written out TODO get rid of
                  // this
  int iring;      // Ring index of the current ring
  // Variables for the hydrogen-bonded 'second' basal ring
  std::vector<int>
      basal2; // Elements are bonded to each element of basal1 in order
  std::vector<int> basal1; // 'First' basal ring
  std::vector<int>
      unOrderedBasal2; // Unordered basal2 ring, the ith element is not
                       // necessarily bonded to the ith element of basal1
  int iatom, jatom;    // Atom numbers (starting from 1), not C++ indices; saved
                       // inside rings
  int natom;           // Neighbour list ID for iatom
  int findAtom;        // Atom number of atomID to find in neighbour list
  bool atomFound; // bool to check if an atomID has been found in the neighbour
                  // list or not
  // Variables for ordering basal2
  // After finding the nearest neighbours, elements which are not nearest
  // neighbours are assigned a value of zero.
  int needle;      // First non-zero element of basal2 after getting the nearest
                   // neighbours
  int startNeedle; // Index of basal2; the first non-zero element of basal2
  int startHayStack; // Index of unOrderedBasal2, corresponding to the first
                     // non-zero element of basal2
  bool isClock;      // Original order of unOrderedBasal2
  int nextElement;   // Next element in unOrderedBasal2 after startHayStack

  // ----------------
  // Create output dir if it doesn't exist already
  const char *path = "../output/cages"; // relative to the build directory
  fs::path dir(path);
  // if (fs::create_directory(dir)) {
  //   std::cerr << "Output Cage directory created\n";
  // }
  // ----------------
  // Subdirectory

  const fs::path path1("../output/cages/hexBasalRings");
  // fs::create_directories(path1);

  // Write output to file inside the output directory
  std::string fileOutNameFull = "../output/cages/hexBasalRings/" + filename;
  outputFile.open(fileOutNameFull);

  // Format:
  // Coordinate IDs (starting from 1), ordered according to the input XYZ file
  // The first line is a comment line
  outputFile << "# Particle IDs in the two basal rings\n";

  // ---------------
  // Find the nearest neighbours of basal1 elements in basal2
  basal1 = rings[currentCage[0]];          // First basal ring
  unOrderedBasal2 = rings[currentCage[1]]; // Unordered second basal ring

  for (int i = 0; i < basal1.size(); i++) {
    iatom = basal1[i]; // This is the atom particle ID, not the C++ index

    // Search through unOrderedBasal2 for an element in the neighbourlist of
    // iatom
    for (int k = 0; k < unOrderedBasal2.size(); k++) {
      findAtom =
          unOrderedBasal2[k]; // Atom ID to find in the neighbour list of iatom

      atomFound = false; // init
      jatom = 0;

      // Search through the neighbour list for findAtom
      for (int n = 1; n < nList[iatom - 1].size(); n++) {
        natom = nList[iatom - 1][n]; // Atom ID

        if (findAtom == natom) {
          atomFound = true;
          break;
        } // Check
      }   // Loop through nList

      if (atomFound) {
        jatom = natom;
        break;
      } // atom has been found
    }   // end of loop through all atomIDs in unOrderedBasal2
    basal2.push_back(jatom);
  } // end of loop through all the atomIDs in basal1

  // ---------------------------------------------------
  // Get particles which are not nearest neighbours
  // 'Alternately' ordered particles

  // ---------------
  // Init
  isClock = false;

  // ---------------
  // Get the first non-zero index {needle}
  for (int i = 0; i < 2; i++) {
    if (basal2[i] != 0) {
      needle = basal2[i]; // Set the needle to the first non-zero element
      startNeedle = i;    // Index of basal2
      break;              // Break out of the loop
    }                     // end of checking for non-zero index
  }                       // end of getting the first non-zero index

  // Find the index matching needle in unOrderedBasal2
  for (int i = 0; i < unOrderedBasal2.size(); i++) {
    if (unOrderedBasal2[i] == needle) {
      startHayStack = i; // Index at which needle has been found
    }                    // end of check for needle
  }                      // end of search for needle in unOrderedBasal2

  // ---------------
  // Check for 'clockwise' order
  // Check 'next' element
  nextElement = startHayStack + 2;
  if (nextElement >= ringSize) {
    nextElement -= ringSize;
  }

  // Init (indices of basal2 and unOrderedBasal2 respectively)
  iatom = 0;
  jatom = 0;

  if (basal2[startNeedle + 2] == unOrderedBasal2[nextElement]) {
    isClock = true;
    // Fill the three elements in increasing order
    for (int i = 1; i < ringSize; i += 2) {
      iatom = startNeedle + i;
      jatom = startHayStack + i;
      // Make sure the indices are not larger than 6
      if (iatom >= ringSize) {
        iatom -= ringSize;
      }
      if (jatom >= ringSize) {
        jatom -= ringSize;
      }

      basal2[iatom] = unOrderedBasal2[jatom];
    } // end of filling next two alternate elements
  }   // check to see if clockwise order is correct

  // ---------------
  // Check for 'anticlockwise' order
  // Check 'next' element
  if (!isClock) {
    iatom = 0;
    jatom = 0; // init

    // First element
    iatom = startNeedle + 2;
    jatom = startHayStack - 2;
    if (jatom < 0) {
      jatom += ringSize;
    }

    if (basal2[iatom] == unOrderedBasal2[jatom]) {
      // Fill the three elements in increasing order
      for (int i = 1; i < ringSize; i += 2) {
        iatom = startNeedle + i;
        jatom = startHayStack - i;
        // Make sure the indices are not larger than 6
        if (iatom > ringSize) {
          iatom -= ringSize;
        }
        if (jatom < 0) {
          jatom += ringSize;
        }

        basal2[iatom] = unOrderedBasal2[jatom];
      } // end of filling next two alternate elements
    }   // check to see if anticlockwise order is correct
    else {
      std::cerr << "Something is wrong with your HCs.\n";
      return 1;
    } // exit with error
  }   // End of check for anticlockwise stuff

  // ---------------------------------------------------
  // Print out the ordered rings
  // For hexagonal cages:
  // Only print out basal1 and basal2

  // BASAL1
  for (int i = 0; i < basal1.size(); i++) {
    outputFile << basal1[i] << " ";
  } // end of loop through basal1
  outputFile << "\n";

  // BASAL2
  for (int i = 0; i < basal2.size(); i++) {
    outputFile << basal2[i] << " ";
  } // end of loop through basal2

  // Close the output file
  outputFile.close();

  return 0;
}

/**
 * @details Function for printing out the basal rings of the prism blocks
 * described by the number prismNum. Uses Boost!
 */
int sout::writeBasalRingsPrism(
    std::vector<int> *basal1, std::vector<int> *basal2, int prismNum,
    std::vector<std::vector<int>> nList,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    bool isDeformed) {
  std::ofstream outputFile;
  std::string number = std::to_string(prismNum);
  std::string filename = "basalRings" + number + ".dat";
  int ringSize =
      (*basal1).size(); // Size of the ring; each ring contains n elements
  int nBonds;           // Number of bonds between the two deformed prisms
  int l_k, m_k;         // Atom ID in basal1 and basal2 respectively
  int iatom,
      jatom; // Index not ID of basal1 and basal2 respectively which match
  int currentIatom, currentJatom; // Index not ID for matchedBasal1 and
                                  // matchedBasal2 respectively
  bool isNeighbour; // Basal rings should have at least one nearest neighbour
  bool isClock;     // The order is in the original 'direction' of basal2
  bool isAntiClock; // The order must be reversed
  std::vector<int> matchedBasal1, matchedBasal2; // Vectors with revised order

  // ----------------
  // Path for deformed prisms
  if (isDeformed) {
    // Create output dir if it doesn't exist already
    const char *path = "../output/deformed"; // relative to the build directory
    fs::path dir(path);
    // if (fs::create_directory(dir)) {
    //   std::cerr << "Output Cage directory created\n";
    // }
    // ----------------
    // Subdirectory

    const fs::path path1("../output/deformed/basalRings");
    // fs::create_directories(path1);

    // Write output to file inside the output directory
    std::string fileOutNameFull = "../output/deformed/basalRings/" + filename;
    outputFile.open(fileOutNameFull);
  } else {
    // Create output dir if it doesn't exist already
    const char *path = "../output/perfect"; // relative to the build directory
    fs::path dir(path);
    // if (fs::create_directory(dir)) {
    //   std::cerr << "Output Cage directory created\n";
    // }
    // ----------------
    // Subdirectory

    const fs::path path1("../output/perfect/basalRings");
    // fs::create_directories(path1);

    // Write output to file inside the output directory
    std::string fileOutNameFull = "../output/perfect/basalRings/" + filename;
    outputFile.open(fileOutNameFull);
  } // end of creating file paths

  // Format:
  // Coordinate IDs (starting from 1), ordered according to the input XYZ file
  // The first line is a comment line
  outputFile << "# Particle IDs in the two basal rings\n";

  // ---------------
  // Find the nearest neighbours of basal1 elements in basal2

  nBonds = 0;
  isNeighbour = false;
  // Loop through every element of basal1
  for (int l = 0; l < ringSize; l++) {
    l_k = (*basal1)[l]; // This is the atom particle ID, not the C++ index

    // Search for the nearest neighbour of l_k in basal2
    // Loop through basal2 elements
    for (int m = 0; m < ringSize; m++) {
      m_k = (*basal2)[m]; // Atom ID to find in the neighbour list of iatom

      // Find m_k inside l_k neighbour list
      auto it = std::find(nList[l_k].begin() + 1, nList[l_k].end(), m_k);

      // If the element has been found, for l1
      if (it != nList[l_k].end()) {
        isNeighbour = true;
        iatom = l; // index of basal1
        jatom = m; // index of basal2
        break;
      } // found element

    } // end of loop through all atomIDs in basal2

    if (isNeighbour) {
      break;
    } // nearest neighbour found
  }   // end of loop through all the atomIDs in basal1

  if (!isNeighbour) {
    std::cerr << "Something is wrong with your deformed prism.\n";
    return 1;
  }
  // ---------------------------------------------------
  // Find out if the order of basal2 is 'clockwise' or 'anticlockwise'
  isClock = false; // init
  isAntiClock = false;

  // atom index (not ID)
  int tempJfor, tempJback;

  tempJfor = jatom + 1;
  tempJback = jatom - 1;

  if (jatom == ringSize - 1) {
    tempJfor = 0;
    tempJback = ringSize - 2;
  }
  if (jatom == 0) {
    tempJfor = 1;
    tempJback = ringSize - 1;
  }

  int forwardJ = (*basal2)[tempJfor];
  int backwardJ = (*basal2)[tempJback];
  int currentI = (*basal1)[iatom];

  // Check clockwise
  double distClock = gen::periodicDist(yCloud, currentI, forwardJ);
  double distAntiClock = gen::periodicDist(yCloud, currentI, backwardJ);

  // Clockwise
  if (distClock < distAntiClock) {
    isClock = true;
  } // end of clockwise check
  // Anti-clockwise
  if (distAntiClock < distClock) {
    isAntiClock = true;
  } // end of anti-clockwise check
  // Some error
  if (isClock == false && isAntiClock == false) {
    // std::cerr << "The points are equidistant.\n";
    return 1;
  } // end of error handling
  // ---------------------------------------------------
  // Get the order of basal1 and basal2
  for (int i = 0; i < ringSize; i++) {
    currentIatom = iatom + i;
    if (currentIatom >= ringSize) {
      currentIatom -= ringSize;
    } // end of basal1 element wrap-around

    // In clockwise order
    if (isClock) {
      currentJatom = jatom + i;
      if (currentJatom >= ringSize) {
        currentJatom -= ringSize;
      } // wrap around
    }   // end of clockwise update
    else {
      currentJatom = jatom - i;
      if (currentJatom < 0) {
        currentJatom += ringSize;
      } // wrap around
    }   // end of anti-clockwise update

    // Add to matchedBasal1 and matchedBasal2 now
    matchedBasal1.push_back((*basal1)[currentIatom]);
    matchedBasal2.push_back((*basal2)[currentJatom]);
  }
  // ---------------------------------------------------
  // Print out the ordered rings
  // For hexagonal cages:
  // Only print out basal1 and basal2

  // BASAL1
  for (int i = 0; i < matchedBasal1.size(); i++) {
    outputFile << matchedBasal1[i] << " ";
  } // end of loop through basal1
  outputFile << "\n";

  // BASAL2
  for (int i = 0; i < matchedBasal2.size(); i++) {
    outputFile << matchedBasal2[i] << " ";
  } // end of loop through basal2

  // Close the output file
  outputFile.close();

  return 0;
}

/**
 * @details Function for printing out cluster statistics
 */
int sout::writeClusterStats(std::string path, int currentFrame,
                            int largestCluster, int numOfClusters,
                            int smallestCluster, double avgClusterSize,
                            int firstFrame) {
  std::ofstream outputFile;
  // ----------------
  // Make the output directory if it doesn't exist
  sout::makePath(path);
  // ----------------
  // Write output to file inside the output directory
  outputFile.open(path + "clusterStats.dat",
                  std::ios_base::app | std::ios_base::out);

  // Format:
  // Comment line
  // 1 3 0 0 4 35 40 ....

  // ----------------
  // Comment line for the first frame
  if (currentFrame == firstFrame) {
    outputFile << "Frame largestCluster numOfClusters smallestCluster "
                  "avgClusterSize\n";
  }
  // ----------------

  outputFile << currentFrame << " " << largestCluster << " " << numOfClusters
             << " " << smallestCluster << " " << avgClusterSize << "\n";

  outputFile.close();

  return 0;
}

/**
 * @details Function for printing out the molecule IDs present in the slice.
 The format should be compatible with the group command in LAMMPS 
 */
int sout::writeMoleculeIDsInSlice(std::string path,
                            molSys::PointCloud<molSys::Point<double>, double> *yCloud) {
  std::ofstream outputFile;
  std::string filename =
      "molID-" + std::to_string(yCloud->currentFrame) + ".dat";
  std::vector<int> idVec; // Vector which will contain all molecule IDs present in the slice
  // Eventually we will sort this and only keep unique molecule IDs. 
  int prevElem, currentElem; // previous and current mol ID in the sequence
  int lastPrintedElemSeries; // Last element printed 
  // ----------------
  // Make the output directory if it doesn't exist
  sout::makePath(path);
  sout::makePath(path+"selection");
  sout::makePath(path+"selection/IDtextFiles");
  // ----------------
  // Write output to file inside the output directory
  outputFile.open(path + "selection/IDtextFiles/" + filename);

  // ----------------
  // Find all molecule IDs from yCloud 
  // Loop through all iatom in yCloud
  for (int iatom = 0; iatom < yCloud->nop; iatom++)
  {
    // If iatom is in the slice, add the molecule ID to idVec
    if (yCloud->pts[iatom].inSlice)
    {
      idVec.push_back(yCloud->pts[iatom].molID);
    } // end of adding molecule ID to vector for slice 
  } // end of loop through iatom
  // ----------------
  // Sort in ascending order
  std::sort(idVec.begin(), idVec.end());
  auto it = std::unique(idVec.begin(), idVec.end());
  // Get rid of undefined elements
  idVec.resize(std::distance(idVec.begin(), it));
  // ----------------

  // Format:
  // Comment line
  // 1 2 3 4 6 7

  // ----------------
  // Comment line 
  outputFile << "# Molecule IDs in slice\n";
  outputFile << "# LAMMPS command : group groupName molecule 100:10000 \n";
  // Format
  // In LAMMPS, groups can be assigned by ID using the following command 
  // group groupName molecule 100:10000 
  // ----------------
  // First element 
  outputFile << idVec[0];
  prevElem = idVec[0];
  lastPrintedElemSeries = idVec[0];
  // ----------------

  // Print other molecule IDs to the file
  for (int i=1; i<idVec.size(); i++)
  {
    currentElem = idVec[i]; // current mol ID
    // prevElem is the previous mol ID
    // if the currentElem-prevElem>1
    if (currentElem-prevElem>1 || i==idVec.size()-1)
    {
      if (lastPrintedElemSeries!=prevElem)
      {
        outputFile << ":" << prevElem << " " << currentElem;
        lastPrintedElemSeries = currentElem;
      } // the previous element has not been printed
      else{
        outputFile << " " << currentElem;
        lastPrintedElemSeries = currentElem;
      } // the previous element has already been printed 
      //
    } // print currentElem
    //
    prevElem = currentElem;
  } // end of printing all elements to file 

  outputFile.close();

  return 0;
}

/**
 * @details Function for printing out the molecule IDs present in the slice.
 The format should be compatible with the group command in LAMMPS 
 */
int sout::writeMoleculeIDsExpressionSelectOVITO(std::string path,
                            molSys::PointCloud<molSys::Point<double>, double> *yCloud) {
  std::ofstream outputFile;
  std::string filename =
      "ovito-molIDSelect-" + std::to_string(yCloud->currentFrame) + ".dat";
  std::vector<int> idVec; // Vector which will contain all molecule IDs present in the slice
  // Eventually we will sort this and only keep unique molecule IDs. 
  int prevElem, currentElem; // previous and current mol ID in the sequence
  int lastPrintedElemSeries; // Last element printed 
  // ----------------
  // Make the output directory if it doesn't exist
  sout::makePath(path);
  sout::makePath(path+"selection");
  sout::makePath(path+"selection/IDovitoFiles");
  // ----------------
  // Write output to file inside the output directory
  outputFile.open(path + "selection/IDovitoFiles/" + filename);

  // ----------------
  // Find all molecule IDs from yCloud 
  // Loop through all iatom in yCloud
  for (int iatom = 0; iatom < yCloud->nop; iatom++)
  {
    // If iatom is in the slice, add the molecule ID to idVec
    if (yCloud->pts[iatom].inSlice)
    {
      idVec.push_back(yCloud->pts[iatom].molID);
    } // end of adding molecule ID to vector for slice 
  } // end of loop through iatom
  // ----------------
  // Sort in ascending order
  std::sort(idVec.begin(), idVec.end());
  auto it = std::unique(idVec.begin(), idVec.end());
  // Get rid of undefined elements
  idVec.resize(std::distance(idVec.begin(), it));
  // ----------------

  // Format:
  // Comment line
  // 1 2 3 4 6 7

  // ----------------
  // Comment line 
  outputFile << "# Molecule IDs in slice\n";
  outputFile << "# OVITO Expression select command \n";
  // ----------------

  // Print other molecule IDs to the file
  for (int i=0; i<idVec.size()-1; i++)
  {
    currentElem = idVec[i]; // current mol ID

    outputFile << "MoleculeIdentifier == " << currentElem << " || ";

  } // end of printing all elements to file except the last one 

  // Print the last element
  outputFile << "MoleculeIdentifier == " << idVec.back();

  outputFile.close();

  return 0;
}

/**
 * @details Function for printing out prism info, when there is no  volume slice
 */
int sout::writePrismNum(std::string path, std::vector<int> nPrisms,
                        std::vector<int> nDefPrisms,
                        std::vector<double> heightPercent, int maxDepth,
                        int currentFrame, int firstFrame) {
  std::ofstream outputFile;
  int totalPrisms; // Number of total prisms
  // ----------------
  // Make the output directory if it doesn't exist
  sout::makePath(path);
  std::string outputDirName = path + "topoINT";
  sout::makePath(outputDirName);
  // ----------------
  // Write output to file inside the output directory
  outputFile.open(path + "topoINT/nPrisms.dat",
                  std::ios_base::app | std::ios_base::out);

  // ----------------
  // Write the comment line if the first frame is being written out
  if (currentFrame == firstFrame) {
    outputFile << "Frame RingSize Num_of_prisms Height% RingSize ... Height\n";
  }
  // ----------------
  // Format:
  // Frame RingSize Num_of_prisms Height% RingSize ... Height%
  // 1 3 0 0 4 35 40 ....

  outputFile << currentFrame << " ";

  for (int ringSize = 3; ringSize <= maxDepth; ringSize++) {
    totalPrisms = nPrisms[ringSize - 3] + nDefPrisms[ringSize - 3];
    // Write out
    outputFile << ringSize << " " << totalPrisms << " "
               << nDefPrisms[ringSize - 3] << " " << heightPercent[ringSize - 3]
               << " ";
  }

  outputFile << "\n";

  outputFile.close();

  return 0;
}

/**
 * @details Function for printing out ring info, for a monolayer
 */
int sout::writeRingNum(std::string path, int currentFrame,
                       std::vector<int> nRings,
                       std::vector<double> coverageAreaXY,
                       std::vector<double> coverageAreaXZ,
                       std::vector<double> coverageAreaYZ, int maxDepth,
                       int firstFrame) {
  std::ofstream outputFileXY;
  std::ofstream outputFileXZ;
  std::ofstream outputFileYZ;
  // ----------------
  // Make the output directory if it doesn't exist
  sout::makePath(path);
  std::string outputDirName = path + "topoMonolayer";
  sout::makePath(outputDirName);
  // ----------------
  // Coverage Area of XY
  // Write output to file inside the output directory
  outputFileXY.open(path + "topoMonolayer/coverageAreaXY.dat",
                    std::ios_base::app | std::ios_base::out);

  // Format:
  // Comment line
  // 1 3 0 0 4 35 40 ....

  // ----------------
  // Add comment for the first frame
  if (currentFrame == firstFrame) {
    outputFileXY << "Frame RingSize Num_of_rings CoverageAreaXY% RingSize ... "
                    "CoverageAreaXY%\n";
  }
  // ----------------

  outputFileXY << currentFrame << " ";

  for (int ringSize = 3; ringSize <= maxDepth; ringSize++) {
    outputFileXY << ringSize << " " << nRings[ringSize - 3] << " "
                 << coverageAreaXY[ringSize - 3] << " ";
  }

  outputFileXY << "\n";

  outputFileXY.close();
  // ----------------
  // Coverage Area of XZ
  // Write output to file inside the output directory
  outputFileXZ.open(path + "topoMonolayer/coverageAreaXZ.dat",
                    std::ios_base::app | std::ios_base::out);

  // ----------------
  // Add comment for the first frame
  if (currentFrame == firstFrame) {
    outputFileXZ << "Frame RingSize Num_of_rings CoverageAreaXZ% RingSize ... "
                    "CoverageAreaXZ%\n";
  }
  // ----------------

  // Format:
  // Frame RingSize Num_of_prisms Height% RingSize ... Height%
  // 1 3 0 0 4 35 40 ....

  outputFileXZ << currentFrame << " ";

  for (int ringSize = 3; ringSize <= maxDepth; ringSize++) {
    outputFileXZ << ringSize << " " << nRings[ringSize - 3] << " "
                 << coverageAreaXZ[ringSize - 3] << " ";
  }

  outputFileXZ << "\n";

  outputFileXZ.close();
  // ----------------
  // Coverage Area of YZ
  // Write output to file inside the output directory
  outputFileYZ.open(path + "topoMonolayer/coverageAreaYZ.dat",
                    std::ios_base::app | std::ios_base::out);

  // ----------------
  // Add comment for the first frame
  if (currentFrame == firstFrame) {
    outputFileYZ << "Frame RingSize Num_of_rings CoverageAreaYZ% RingSize ... "
                    "CoverageAreaYZ%\n";
  }
  // ----------------

  // Format:
  // Frame RingSize Num_of_prisms Height% RingSize ... Height%
  // 1 3 0 0 4 35 40 ....

  outputFileYZ << currentFrame << " ";

  for (int ringSize = 3; ringSize <= maxDepth; ringSize++) {
    outputFileYZ << ringSize << " " << nRings[ringSize - 3] << " "
                 << coverageAreaYZ[ringSize - 3] << " ";
  }

  outputFileYZ << "\n";

  outputFileYZ.close();

  return 0;
}

/**
 * @details Function for printing out ring info, for a bulk system
 */
int sout::writeRingNumBulk(std::string path, int currentFrame,
                       std::vector<int> nRings,
                       int maxDepth,
                       int firstFrame) {
  std::ofstream outputFile;
  // ----------------
  // Make the output directory if it doesn't exist
  sout::makePath(path);
  std::string outputDirName = path + "bulkTopo";
  sout::makePath(outputDirName);
  // ----------------
  // Ring output file
  // Write output to file inside the output directory
  outputFile.open(path + "bulkTopo/num_rings.dat",
                    std::ios_base::app | std::ios_base::out);

  // Format:
  // Comment line
  // 1 3 0 4 35 ....

  // ----------------
  // Add comment for the first frame
  if (currentFrame == firstFrame) {
    outputFile << "Frame RingSize Num_of_rings RingSize Num_of_rings...\n";
  }
  // ----------------

  outputFile << currentFrame << " ";

  for (int ringSize = 3; ringSize <= maxDepth; ringSize++) {
    outputFile << ringSize << " " << nRings[ringSize - 3] << " ";
  }

  outputFile << "\n";

  outputFile.close();

  return 0;
}

/**
 * @details Function for printing out the RDF to a file, given that the file
 * already exists and given the filename.
 */
int sout::printRDF(std::string fileName, std::vector<double> *rdfValues,
                   double binwidth, int nbin) {
  //
  std::ofstream outputFile; // For the output file
  double r;                 // Distance for the current bin

  // Append to the file
  outputFile.open(fileName, std::ios_base::app | std::ios_base::out);

  // Loop through all the bins
  for (int ibin = 0; ibin < nbin; ibin++) {
    //
    r = binwidth * (ibin + 0.5); // Current distance for ibin
    outputFile << r << " " << (*rdfValues)[ibin] << "\n";
  } // end of loop through all bins

  outputFile.close();

  return 0;
}

/**
 * @details Function for out the number of DDCs, HCs, mixed rings, basal and
 * prismatic rings.
 */
int sout::writeTopoBulkData(std::string path, int currentFrame, int numHC,
                            int numDDC, int mixedRings, int basalRings,
                            int prismaticRings, int firstFrame) {
  //
  std::ofstream outputFile;
  // ----------------
  // Make the output directory if it doesn't exist
  sout::makePath(path);
  std::string outputDirName = path + "bulkTopo";
  sout::makePath(outputDirName);
  // ----------------
  // Write output to file inside the output directory
  outputFile.open(path + "bulkTopo/cageData.dat",
                  std::ios_base::app | std::ios_base::out);

  // Format:
  // Frame RingSize Num_of_prisms Height% RingSize ... Height%
  // 1 3 0 0 4 35 40 ....
  // -------------------
  // If first frame then write the comment line
  if (currentFrame == firstFrame) {
    outputFile << "Frame HCnumber DDCnumber MixedRingNumber PrismaticRings "
                  "basalRings\n";
  }
  // -------------------
  outputFile << currentFrame << " " << numHC << " " << numDDC << " "
             << mixedRings << " " << prismaticRings << " " << basalRings
             << "\n";

  outputFile.close();
  return 0;
} // end of function

/**
 * @details Prints out a LAMMPS dump file for all the cages in bulk ice, for
 * every frame, printing the RMSD per atom as well
 */
int sout::writeLAMMPSdumpCages(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<double> rmsdPerAtom, std::vector<int> atomTypes,
    std::string path, int firstFrame) {
  std::ofstream outputFile;
  int iatom; // Index, not atom ID
  std::string filename =
      "dump-" + std::to_string(yCloud->currentFrame) + ".lammpstrj";
  // ----------------
  // Make the output directory if it doesn't exist
  std::string outputDirName = path + "bulkTopo/dumpFiles";
  sout::makePath(outputDirName);
  // ----------------
  // Write out information about the data types
  if (yCloud->currentFrame == firstFrame) {
    outputFile.open(path + "bulkTopo/typeInfo.dat");
    outputFile << "Atom types in the dump files are:\n";
    outputFile << " Type 0 (dummy) = unidentified phase\n";
    outputFile << " Type 1 (hc) = atom belonging to a Hexagonal Cage.\n";
    outputFile << " Type 2 (ddc) = atom belonging to a Double-Diamond Cage\n";
    outputFile << " Type 3 (mixed) = atom belonging to a mixed ring shared by "
                  "a DDC and HC\n";
    outputFile
        << " Type 4 (pnc) = atom belonging to a pair of pentagonal rings\n";
    outputFile << " Type 5 (mixed2) = atom belonging to a pentagonal "
                  "nanochannel, shared by DDCs/HCs\n";
    outputFile.close();
  } // end of writing out information
  // ----------------
  // Write output to file inside the output directory
  outputFile.open(path + "bulkTopo/dumpFiles/" + filename);
  // ----------------------------------------------------
  // Header Format

  // ITEM: TIMESTEP
  // 0
  // ITEM: NUMBER OF ATOMS
  // 500
  // ITEM: BOX BOUNDS pp pp pp
  // -9.0400100000000005e-01 1.7170999999999999e+01
  // -9.0400100000000005e-01 1.7170999999999999e+01
  // -9.0400100000000005e-01 1.7170999999999999e+01
  // ITEM: ATOMS id mol type x y z rmsd

  // -----------------
  // -------
  // Write the header
  // ITEM: TIMESTEP
  outputFile << "ITEM: TIMESTEP\n";
  // Write out frame number
  outputFile << yCloud->currentFrame << "\n";
  // ITEM: NUMBER OF ATOMS
  outputFile << "ITEM: NUMBER OF ATOMS\n";
  // Number of atoms
  outputFile << yCloud->pts.size() << "\n";
  // ITEM: BOX BOUNDS pp pp pp
  outputFile << "ITEM: BOX BOUNDS pp pp pp\n";
  // Box lengths
  outputFile << yCloud->boxLow[0] << " " << yCloud->boxLow[0] + yCloud->box[0]
             << "\n";
  outputFile << yCloud->boxLow[1] << " " << yCloud->boxLow[1] + yCloud->box[1]
             << "\n";
  outputFile << yCloud->boxLow[2] << " " << yCloud->boxLow[2] + yCloud->box[2]
             << "\n";
  // ITEM: ATOMS id mol type x y z rmsd
  outputFile << "ITEM: ATOMS id mol type x y z rmsd\n";
  // -------
  // Write out the atom coordinates
  // Format
  // ITEM: ATOMS id mol type x y z rmsd
  //
  // Loop through atoms
  for (int i = 0; i < yCloud->pts.size(); i++) {
    iatom =
        yCloud->pts[i].atomID; // The actual ID can be different from the index
    // Write out coordinates
    outputFile << iatom << " " << yCloud->pts[i].molID << " " << atomTypes[i]
               << " " << yCloud->pts[i].x << " " << yCloud->pts[i].y << " "
               << yCloud->pts[i].z << " " << rmsdPerAtom[i] << "\n";

  } // end of loop through all atoms in pointCloud
  // -----------------------------------------------------
  outputFile.close(); // Close the file
  return 0;
} // end of function

/**
 * @details Prints out a LAMMPS dump file for all the prisms, for every frame,
 * printing the RMSD per atom as well
 */
int sout::writeLAMMPSdumpINT(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<double> rmsdPerAtom, std::vector<int> atomTypes, int maxDepth,
    std::string path) {
  //
  std::ofstream outputFile;
  int iatom; // Index, not atom ID
  std::string filename =
      "dump-" + std::to_string(yCloud->currentFrame) + ".lammpstrj";
  // ----------------
  // Make the output directory if it doesn't exist
  std::string outputDirName = path + "topoINT/dumpFiles";
  sout::makePath(outputDirName);
  // ----------------
  // Write output to file inside the output directory
  outputFile.open(path + "topoINT/dumpFiles/" + filename);
  // ----------------------------------------------------
  // Header Format

  // ITEM: TIMESTEP
  // 0
  // ITEM: NUMBER OF ATOMS
  // 500
  // ITEM: BOX BOUNDS pp pp pp
  // -9.0400100000000005e-01 1.7170999999999999e+01
  // -9.0400100000000005e-01 1.7170999999999999e+01
  // -9.0400100000000005e-01 1.7170999999999999e+01
  // ITEM: ATOMS id mol type x y z rmsd

  // -----------------
  // -------
  // Write the header
  // ITEM: TIMESTEP
  outputFile << "ITEM: TIMESTEP\n";
  // Write out frame number
  outputFile << yCloud->currentFrame << "\n";
  // ITEM: NUMBER OF ATOMS
  outputFile << "ITEM: NUMBER OF ATOMS\n";
  // Number of atoms
  outputFile << yCloud->pts.size() << "\n";
  // ITEM: BOX BOUNDS pp pp pp
  outputFile << "ITEM: BOX BOUNDS pp pp pp\n";
  // Box lengths
  outputFile << yCloud->boxLow[0] << " " << yCloud->boxLow[0] + yCloud->box[0]
             << "\n";
  outputFile << yCloud->boxLow[1] << " " << yCloud->boxLow[1] + yCloud->box[1]
             << "\n";
  outputFile << yCloud->boxLow[2] << " " << yCloud->boxLow[2] + yCloud->box[2]
             << "\n";
  // ITEM: ATOMS id mol type x y z rmsd
  outputFile << "ITEM: ATOMS id mol type x y z rmsd\n";
  // -------
  // Write out the atom coordinates
  // Format
  // ITEM: ATOMS id mol type x y z rmsd
  //
  // Loop through atoms
  for (int i = 0; i < yCloud->pts.size(); i++) {
    iatom =
        yCloud->pts[i].atomID; // The actual ID can be different from the index
    // Write out coordinates
    outputFile << iatom << " " << yCloud->pts[i].molID << " " << atomTypes[i]
               << " " << yCloud->pts[i].x << " " << yCloud->pts[i].y << " "
               << yCloud->pts[i].z << " " << rmsdPerAtom[i] << "\n";

  } // end of loop through all atoms in pointCloud
  // -----------------------------------------------------
  return 0;
} // end of function

/**
 * @details Prints out a LAMMPS dump file for all atoms, for every frame,
 * printing the inSlice attribute for a user-defined slice in a separate column 
 */
int sout::writeLAMMPSdumpSlice(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::string path) {
  //
  std::ofstream outputFile;
  int iatom; // Index, not atom ID
  std::string filename =
      "dump-" + std::to_string(yCloud->currentFrame) + ".lammpstrj";
  // ----------------
  // Make the output directory if it doesn't exist
  sout::makePath(path+"selection");
  std::string outputDirName = path + "selection/dumpFiles";
  sout::makePath(outputDirName);
  // ----------------
  // Write output to file inside the output directory
  outputFile.open(path + "selection/dumpFiles/" + filename);
  // ----------------------------------------------------
  // Header Format

  // ITEM: TIMESTEP
  // 0
  // ITEM: NUMBER OF ATOMS
  // 500
  // ITEM: BOX BOUNDS pp pp pp
  // -9.0400100000000005e-01 1.7170999999999999e+01
  // -9.0400100000000005e-01 1.7170999999999999e+01
  // -9.0400100000000005e-01 1.7170999999999999e+01
  // ITEM: ATOMS id mol type x y z rmsd

  // -----------------
  // -------
  // Write the header
  // ITEM: TIMESTEP
  outputFile << "ITEM: TIMESTEP\n";
  // Write out frame number
  outputFile << yCloud->currentFrame << "\n";
  // ITEM: NUMBER OF ATOMS
  outputFile << "ITEM: NUMBER OF ATOMS\n";
  // Number of atoms
  outputFile << yCloud->pts.size() << "\n";
  // ITEM: BOX BOUNDS pp pp pp
  outputFile << "ITEM: BOX BOUNDS pp pp pp\n";
  // Box lengths
  outputFile << yCloud->boxLow[0] << " " << yCloud->boxLow[0] + yCloud->box[0]
             << "\n";
  outputFile << yCloud->boxLow[1] << " " << yCloud->boxLow[1] + yCloud->box[1]
             << "\n";
  outputFile << yCloud->boxLow[2] << " " << yCloud->boxLow[2] + yCloud->box[2]
             << "\n";
  // ITEM: ATOMS id mol type x y z rmsd
  outputFile << "ITEM: ATOMS id mol type x y z inSlice\n";
  // -------
  // Write out the atom coordinates
  // Format
  // ITEM: ATOMS id mol type x y z rmsd
  //
  // Loop through atoms
  for (int i = 0; i < yCloud->pts.size(); i++) {
    iatom =
        yCloud->pts[i].atomID; // The actual ID can be different from the index
    // Write out coordinates
    outputFile << iatom << " " << yCloud->pts[i].molID << " " << yCloud->pts[i].type
               << " " << yCloud->pts[i].x << " " << yCloud->pts[i].y << " "
               << yCloud->pts[i].z << " " << yCloud->pts[i].inSlice << "\n";

  } // end of loop through all atoms in pointCloud
  // -----------------------------------------------------
  return 0;
} // end of function

/**
 * @details Prints out a LAMMPS data file for all the prisms for a single frame,
 * with some default options. Only Oxygen atoms are printed out. Bonds are
 * inferred from the rings vector of vectors
 */
int sout::writeLAMMPSdataAllPrisms(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> nList, std::vector<int> atomTypes,
    int maxDepth, std::string path, bool doShapeMatching) {
  //
  std::ofstream outputFile;
  int iatom; // Index, not atom ID
  int bondTypes = 1;
  // Bond stuff
  std::vector<std::vector<int>> bonds; // Vector of vector, with each row
                                       // containing the atom IDs of each bond
  std::string filename =
      "system-prisms-" + std::to_string(yCloud->currentFrame) + ".data";

  // ---------------
  // Get the bonds
  bonds = bond::populateBonds(nList, yCloud);
  //
  // ----------------
  // Make the output directory if it doesn't exist
  sout::makePath(path);
  std::string outputDirName = path + "topoINT";
  sout::makePath(outputDirName);
  outputDirName = path + "topoINT/dataFiles/";
  sout::makePath(outputDirName);
  // ----------------
  // Write output to file inside the output directory
  outputFile.open(path + "topoINT/dataFiles/" + filename);
  // FORMAT:
  //  Comment Line
  //  4 atoms
  //  4 bonds
  //  0 angles
  //  0 dihedrals
  //  0 impropers
  //  1 atom types
  //  1 bond types
  //  0 angle types
  //  0 dihedral types
  //  0 improper types
  //  -1.124000 52.845002  xlo xhi
  //  0.000000 54.528999  ylo yhi
  //  1.830501 53.087501  zlo zhi

  //  Masses

  //  1 15.999400 # O

  //  Atoms

  // 1 1 1 0 20.239  1.298 6.873 # O
  // 2 1 1 0 0 5.193 6.873 # O
  // 3 1 1 0 2.249 1.298 6.873 # O

  // -------
  // Write the header
  // Write comment line
  outputFile << "Written out by D-SEAMS\n";
  // Write out the number of atoms
  outputFile << yCloud->pts.size() << " "
             << "atoms"
             << "\n";
  // Number of bonds
  outputFile << bonds.size() << " bonds"
             << "\n";
  outputFile << "0 angles\n0 dihedrals\n0 impropers\n";
  // There are maxDepth-2 total types of prisms + dummy
  if (doShapeMatching) {
    outputFile << 2 * maxDepth - 2 << " atom types\n";
  } else {
    outputFile << maxDepth << " atom types\n";
  }
  // Bond types
  outputFile
      << bondTypes
      << " bond types\n0 angle types\n0 dihedral types\n0 improper types\n";
  // Box lengths
  outputFile << yCloud->boxLow[0] << " " << yCloud->boxLow[0] + yCloud->box[0]
             << " xlo xhi\n";
  outputFile << yCloud->boxLow[1] << " " << yCloud->boxLow[1] + yCloud->box[1]
             << " ylo yhi\n";
  outputFile << yCloud->boxLow[2] << " " << yCloud->boxLow[2] + yCloud->box[2]
             << " zlo zhi\n";
  // Masses
  outputFile << "\nMasses\n\n";
  outputFile << "1 15.999400 # dummy\n";
  outputFile << "2 1.0 # mixedRings \n";
  // There are maxDepth-2 other prism types
  for (int ringSize = 3; ringSize <= maxDepth; ringSize++) {
    outputFile << ringSize << " 15.999400 # prism" << ringSize << "\n";
  } // end of writing out perfect atom types
  // Write out the types for the deformed prism blocks
  if (doShapeMatching) {
    for (int ringSize = maxDepth + 1; ringSize <= 2 * maxDepth - 2;
         ringSize++) {
      int p = ringSize - maxDepth + 2;
      outputFile << ringSize << " 15.999400 # deformPrism" << p << "\n";
    } // end of writing out perfect atom types
  }   // Deformed prism types
  // Atoms
  outputFile << "\nAtoms\n\n";
  // -------
  // Write out the atom coordinates
  // Loop through atoms
  for (int i = 0; i < yCloud->pts.size(); i++) {
    iatom =
        yCloud->pts[i].atomID; // The actual ID can be different from the index
    // Write out coordinates
    // atomID molecule-tag atom-type q x y z
    outputFile << iatom << " " << yCloud->pts[i].molID << " " << atomTypes[i]
               << " 0 " << yCloud->pts[i].x << " " << yCloud->pts[i].y << " "
               << yCloud->pts[i].z << "\n";

  } // end of loop through all atoms in pointCloud

  // Print the bonds now!
  outputFile << "\nBonds\n\n";
  // Loop through all bonds
  for (int ibond = 0; ibond < bonds.size(); ibond++) {
    //
    outputFile << ibond + 1 << " 1 " << bonds[ibond][0] << " "
               << bonds[ibond][1] << "\n";
  } // end of for loop for bonds

  outputFile.close();
  // Once the datafile has been printed, exit
  return 0;
}

/**
 * @details Prints out a LAMMPS data file for all the rings for a single frame,
 * with some default options. Only Oxygen atoms are printed out. Bonds are
 * inferred from the rings vector of vectors
 */
int sout::writeLAMMPSdataAllRings(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> nList, std::vector<int> atomTypes,
    int maxDepth, std::string path, bool isMonolayer) {
  //
  std::ofstream outputFile;
  int iatom; // Index, not atom ID
  int bondTypes = 1;
  // Bond stuff
  std::vector<std::vector<int>> bonds; // Vector of vector, with each row
                                       // containing the atom IDs of each bond
  std::string filename =
      "system-rings-" + std::to_string(yCloud->currentFrame) + ".data";
  std::string pathName, pathFolder;

  // ---------------
  // Get the bonds
  bonds = bond::populateBonds(nList, yCloud);
  //
  // ----------------
  // Make the output directory if it doesn't exist
  if (isMonolayer)
  {
    pathFolder = "topoMonolayer";
    pathName = "topoMonolayer/dataFiles/";
  } else{
    pathFolder = "bulkTopo";
    pathName = "bulkTopo/dataFiles/";
  }
  
  sout::makePath(path);
  std::string outputDirName = path + pathFolder;
  sout::makePath(outputDirName);
  outputDirName = path + pathName;
  sout::makePath(outputDirName);

  // Write output to file inside the output directory
  outputFile.open(path + pathName + filename);

  // FORMAT:
  //  Comment Line
  //  4 atoms
  //  4 bonds
  //  0 angles
  //  0 dihedrals
  //  0 impropers
  //  1 atom types
  //  1 bond types
  //  0 angle types
  //  0 dihedral types
  //  0 improper types
  //  -1.124000 52.845002  xlo xhi
  //  0.000000 54.528999  ylo yhi
  //  1.830501 53.087501  zlo zhi

  //  Masses

  //  1 15.999400 # O

  //  Atoms

  // 1 1 1 0 20.239  1.298 6.873 # O
  // 2 1 1 0 0 5.193 6.873 # O
  // 3 1 1 0 2.249 1.298 6.873 # O

  // -------
  // Write the header
  // Write comment line
  outputFile << "Written out by D-SEAMS\n";
  // Write out the number of atoms
  outputFile << yCloud->pts.size() << " "
             << "atoms"
             << "\n";
  // Number of bonds
  outputFile << bonds.size() << " bonds"
             << "\n";
  outputFile << "0 angles\n0 dihedrals\n0 impropers\n";
  // There are maxDepth-2 total types of prisms + dummy
  outputFile << maxDepth << " atom types\n";
  // Bond types
  outputFile
      << bondTypes
      << " bond types\n0 angle types\n0 dihedral types\n0 improper types\n";
  // Box lengths
  outputFile << yCloud->boxLow[0] << " " << yCloud->boxLow[0] + yCloud->box[0]
             << " xlo xhi\n";
  outputFile << yCloud->boxLow[1] << " " << yCloud->boxLow[1] + yCloud->box[1]
             << " ylo yhi\n";
  outputFile << yCloud->boxLow[2] << " " << yCloud->boxLow[2] + yCloud->box[2]
             << " zlo zhi\n";
  // Masses
  outputFile << "\nMasses\n\n";
  outputFile << "1 15.999400 # dummy\n";
  outputFile << "2 1.0 # \n";
  // There are maxDepth-2 other prism types
  for (int ringSize = 3; ringSize <= maxDepth; ringSize++) {
    outputFile << ringSize << " 15.999400 # ring" << ringSize << "\n";
  } // end of writing out atom types
  // Atoms
  outputFile << "\nAtoms\n\n";
  // -------
  // Write out the atom coordinates
  // Loop through atoms
  for (int i = 0; i < yCloud->pts.size(); i++) {
    iatom =
        yCloud->pts[i].atomID; // The actual ID can be different from the index
    // Write out coordinates
    // atomID molecule-tag atom-type q x y z
    outputFile << iatom << " " << yCloud->pts[i].molID << " " << atomTypes[i]
               << " 0 " << yCloud->pts[i].x << " " << yCloud->pts[i].y << " "
               << yCloud->pts[i].z << "\n";

  } // end of loop through all atoms in pointCloud

  // Print the bonds now!
  outputFile << "\nBonds\n\n";
  // Loop through all bonds
  for (int ibond = 0; ibond < bonds.size(); ibond++) {
    //
    outputFile << ibond + 1 << " 1 " << bonds[ibond][0] << " "
               << bonds[ibond][1] << "\n";
  } // end of for loop for bonds

  outputFile.close();
  // Once the datafile has been printed, exit
  return 0;
}

/**
 * @details Prints out a LAMMPS data file for the prisms, with some default
 * options. Only Oxygen atoms are printed out
 */
int sout::writeLAMMPSdataPrisms(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> rings, bool useBondFile, std::string bondFile,
    std::vector<int> listPrism, std::vector<std::vector<int>> nList,
    std::string filename) {
  std::ofstream outputFile;
  std::vector<int> atoms;         // Holds all atom IDs to print
  int ringSize = rings[0].size(); // Ring size of each ring in rings
  int iatom;                      // Index, not atom ID
  bool padAtoms = false;          // Add extra atoms if the atom IDs are skipped
  int prevAtomID = 0;             // Check for previous atom ID
  int dummyAtoms = 0;             // Number of dummy atoms to fill
  int dummyID;
  int jatom; // Array index is 1 less than the ID (index for dummy atom)
  int bondTypes = 1;
  // Bond stuff
  std::vector<std::vector<int>> bonds; // Vector of vector, with each row
                                       // containing the atom IDs of each bond
  bool atomOne, atomTwo;               // If bond atoms are in the prism or not
  bool isPrismBond;

  // ----------------
  // Return if there are no prisms
  if (listPrism.size() == 0) {
    return 1;
  }

  // ---------------
  // Get the bonds
  if (useBondFile) {
    // Bonds from file
    bonds = sinp::readBonds(bondFile);
  } // get bonds from file
  else {
    bonds = bond::populateBonds(nList, yCloud);
  } // Bonds from rings
  //
  // ----------------
  // Otherwise create file
  // Create output dir if it doesn't exist already
  const char *path = "../output"; // relative to the build directory
  fs::path dir(path);
  // if (fs::create_directory(dir)) {
  //   std::cerr << "Output directory created\n";
  // }
  // ----------------
  // Get all the unique atomIDs of the atoms in the rings of this type
  // Put all atom IDs into one 1-D vector
  size_t total_size{0};
  // Get the total number of atoms (repeated)
  total_size = listPrism.size() * ringSize;
  // Reserve this size inside atoms
  atoms.reserve(total_size);
  // Fill up all these atom IDs
  for (int iring = 0; iring < listPrism.size(); iring++) {
    std::move(rings[listPrism[iring]].begin(), rings[listPrism[iring]].end(),
              std::back_inserter(atoms));
  } // end of loop through all rings in the list

  // Sort the array according to atom ID, which will be needed to get the
  // unique IDs and to remove duplicates
  sort(atoms.begin(), atoms.end());
  // Get the unique atom IDs
  auto ip = std::unique(atoms.begin(), atoms.end());
  // Resize the vector to remove undefined terms
  atoms.resize(std::distance(atoms.begin(), ip));
  // If the number of atoms is less than the total nop, add dummy atoms
  if (atoms.size() != yCloud->nop) {
    padAtoms = true;
    bondTypes = 2;
  }
  // ----------------
  // Write output to file inside the output directory
  outputFile.open("../output/" + filename);
  // FORMAT:
  //  Comment Line
  //  4 atoms
  //  4 bonds
  //  0 angles
  //  0 dihedrals
  //  0 impropers
  //  1 atom types
  //  1 bond types
  //  0 angle types
  //  0 dihedral types
  //  0 improper types
  //  -1.124000 52.845002  xlo xhi
  //  0.000000 54.528999  ylo yhi
  //  1.830501 53.087501  zlo zhi

  //  Masses

  //  1 15.999400 # O

  //  Atoms

  // 1 1 1 0 20.239  1.298 6.873 # O
  // 2 1 1 0 0 5.193 6.873 # O
  // 3 1 1 0 2.249 1.298 6.873 # O

  // -------
  // Write the header
  // Write comment line
  outputFile << "Written out by D-SEAMS\n";
  // Write out the number of atoms
  outputFile << yCloud->pts.size() << " "
             << "atoms"
             << "\n";
  // Number of bonds
  outputFile << bonds.size() << " bonds"
             << "\n";
  outputFile << "0 angles\n0 dihedrals\n0 impropers\n";
  // If padded atoms are required, two atom types will be required
  if (padAtoms) {
    outputFile << "2 atom types\n";
  } else {
    outputFile << "1 atom types\n";
  } // end of atom types
  outputFile
      << bondTypes
      << " bond types\n0 angle types\n0 dihedral types\n0 improper types\n";
  // Box lengths
  outputFile << "0 " << yCloud->box[0] << " xlo xhi\n";
  outputFile << "0 " << yCloud->box[1] << " ylo yhi\n";
  outputFile << "0 " << yCloud->box[2] << " zlo zhi\n";
  // Masses
  outputFile << "\nMasses\n\n";
  outputFile << "1 15.999400 # O\n";
  if (padAtoms) {
    outputFile << "2 1.0 # dummy\n";
  }
  // Atoms
  outputFile << "\nAtoms\n\n";
  // -------
  // Write out the atom coordinates
  // Loop through atoms
  for (int i = 0; i < atoms.size(); i++) {
    iatom = atoms[i] - 1; // The actual index is one less than the ID
    // -----------
    // Pad out
    // Fill in dummy atoms if some have been skipped
    if (atoms[i] != prevAtomID + 1) {
      dummyAtoms = atoms[i] - prevAtomID - 1;
      dummyID = prevAtomID;
      // Loop to write out dummy atoms
      for (int j = 0; j < dummyAtoms; j++) {
        dummyID++;
        jatom = dummyID - 1;
        // 1 molecule-tag atom-type q x y z
        outputFile << dummyID << " " << yCloud->pts[jatom].molID << " 2 0 "
                   << yCloud->pts[jatom].x << " " << yCloud->pts[jatom].y << " "
                   << yCloud->pts[jatom].z << "\n";
      } // end of dummy atom write-out
    }   // end of check for dummy atom printing
    // -----------
    // Write out coordinates
    // 1 molecule-tag atom-type q x y z
    outputFile << atoms[i] << " " << yCloud->pts[iatom].molID << " 1 0 "
               << yCloud->pts[iatom].x << " " << yCloud->pts[iatom].y << " "
               << yCloud->pts[iatom].z << "\n";
    // update the previous atom ID
    prevAtomID = atoms[i];
  } // end of loop through all atoms in atomID

  // Fill in the rest of the dummy atoms
  if (atoms[atoms.size() - 1] != yCloud->nop) {
    //
    for (int id = atoms[atoms.size() - 1] + 1; id <= yCloud->nop; id++) {
      jatom = id - 1;
      outputFile << id << " " << yCloud->pts[jatom].molID << " 2 0 "
                 << yCloud->pts[jatom].x << " " << yCloud->pts[jatom].y << " "
                 << yCloud->pts[jatom].z << "\n";
    } // end of printing out dummy atoms
  }

  // Print the bonds now!
  outputFile << "\nBonds\n\n";
  // Loop through all bonds
  for (int ibond = 0; ibond < bonds.size(); ibond++) {
    // Init
    isPrismBond = false;
    atomOne = false;
    atomTwo = false;
    // --------
    // Check if the bond is in the prism or not
    auto it = std::find(atoms.begin() + 1, atoms.end(), bonds[ibond][0]);
    if (it != atoms.end()) {
      atomOne = true;
    } else if (bonds[ibond][0] == atoms[0]) {
      atomOne = true;
    } else if (bonds[ibond][0] == atoms[atoms.size() - 1]) {
      atomOne = true;
    }

    auto it1 = std::find(atoms.begin() + 1, atoms.end(), bonds[ibond][1]);
    if (it1 != atoms.end()) {
      atomTwo = true;
    } else if (bonds[ibond][1] == atoms[0]) {
      atomTwo = true;
    } else if (bonds[ibond][1] == atoms[atoms.size() - 1]) {
      atomTwo = true;
    }

    if (atomOne == false || atomTwo == false) {
      isPrismBond = false;
    } else {
      isPrismBond = true;
    }
    // --------
    if (isPrismBond) {
      // is inside the prism (type 1)
      outputFile << ibond + 1 << " 1 " << bonds[ibond][0] << " "
                 << bonds[ibond][1] << "\n";
    } else {
      // not inside the prism (type 2)
      outputFile << ibond + 1 << " 2 " << bonds[ibond][0] << " "
                 << bonds[ibond][1] << "\n";
    }

  } // end of for loop for bonds

  outputFile.close();
  // Once the datafile has been printed, exit
  return 0;
}

/**
 *  @details Prints out a LAMMPS data file for the either the DDCs or HCs, with
 * some default options. Only Oxygen atoms are printed out
 */
int sout::writeLAMMPSdataCages(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> rings, std::vector<cage::Cage> *cageList,
    cage::cageType type, int numCages, std::string filename) {
  std::ofstream outputFile;
  std::vector<int> atoms;         // Holds all atom IDs to print
  int ringSize = rings[0].size(); // Ring size of each ring in rings
  int iatom;                      // Index, not atom ID
  bool padAtoms = false;          // Add extra atoms if the atom IDs are skipped
  int prevAtomID = 0;             // Check for previous atom ID
  int dummyAtoms = 0;             // Number of dummy atoms to fill
  int dummyID;
  int iring; // Current ring index (vector index)
  int jatom; // Array index is 1 less than the ID (index for dummy atom)
  int bondTypes = 1;
  std::string actualCageType; // The actual name of the cage types
  // Bond stuff
  std::vector<std::vector<int>> bonds; // Vector of vector, with each row
                                       // containing the atom IDs of each bond
  int nRings = 0;                      // Number of rings

  // ----------------
  // Return if there are no cages at all
  if ((*cageList).size() == 0) {
    return 1;
  }

  // Return if there are no cages of the required type
  if (numCages == 0) {
    return 1;
  }
  // ---------------
  // Get the bonds
  bonds = bond::createBondsFromCages(rings, cageList, type, &nRings);
  //
  // ----------------
  // Otherwise create file
  // Create output dir if it doesn't exist already
  const char *path = "../output"; // relative to the build directory
  fs::path dir(path);
  // if (fs::create_directory(dir)) {
  //   std::cerr << "Output directory created\n";
  // }
  // ----------------
  // Get all the unique atomIDs of the atoms in the rings of this type
  // Put all atom IDs into one 1-D vector
  size_t total_size{0};
  // Get the total number of atoms (repeated)
  total_size = nRings * ringSize;
  // Reserve this size inside atoms
  atoms.reserve(total_size);
  // Fill up all these atom IDs
  // Loop through every cage in cageList
  for (int icage = 0; icage < (*cageList).size(); icage++) {
    // Skip if the cage is of a different type
    if ((*cageList)[icage].type != type) {
      continue;
    }
    // Loop through every ring inside Cage
    for (int k = 0; k < (*cageList)[icage].rings.size(); k++) {
      iring = (*cageList)[icage].rings[k]; // Current ring index
      std::move(rings[iring].begin(), rings[iring].end(),
                std::back_inserter(atoms));
    } // end of loop through every ring in icage

  } // end of loop through all cages in cageList

  // Sort the array according to atom ID, which will be needed to get the
  // unique IDs and to remove duplicates
  sort(atoms.begin(), atoms.end());
  // Get the unique atom IDs
  auto ip = std::unique(atoms.begin(), atoms.end());
  // Resize the vector to remove undefined terms
  atoms.resize(std::distance(atoms.begin(), ip));
  // If the number of atoms is less than the total nop, add dummy atoms
  if (atoms.size() != yCloud->nop) {
    padAtoms = true;
    bondTypes = 1;
  }
  // ----------------
  // Write output to file inside the output directory
  outputFile.open("../output/" + filename);
  // FORMAT:
  //  Comment Line
  //  4 atoms
  //  4 bonds
  //  0 angles
  //  0 dihedrals
  //  0 impropers
  //  1 atom types
  //  1 bond types
  //  0 angle types
  //  0 dihedral types
  //  0 improper types
  //  -1.124000 52.845002  xlo xhi
  //  0.000000 54.528999  ylo yhi
  //  1.830501 53.087501  zlo zhi

  //  Masses

  //  1 15.999400 # O

  //  Atoms

  // 1 1 1 0 20.239  1.298 6.873 # O
  // 2 1 1 0 0 5.193 6.873 # O
  // 3 1 1 0 2.249 1.298 6.873 # O

  // -------
  // Write the header
  // Write comment line
  outputFile << "Written out by D-SEAMS\n";
  // Write out the number of atoms
  outputFile << yCloud->pts.size() << " "
             << "atoms"
             << "\n";
  // Number of bonds
  outputFile << bonds.size() << " bonds"
             << "\n";
  outputFile << "0 angles\n0 dihedrals\n0 impropers\n";
  // If padded atoms are required, two atom types will be required
  if (padAtoms) {
    outputFile << "2 atom types\n";
  } else {
    outputFile << "1 atom types\n";
  } // end of atom types
  outputFile
      << bondTypes
      << " bond types\n0 angle types\n0 dihedral types\n0 improper types\n";
  // Box lengths
  outputFile << "0 " << yCloud->box[0] << " xlo xhi\n";
  outputFile << "0 " << yCloud->box[1] << " ylo yhi\n";
  outputFile << "0 " << yCloud->box[2] << " zlo zhi\n";
  // Masses
  outputFile << "\nMasses\n\n";
  // For DDCs and HCs
  if (type == cage::cageType::HexC) {
    actualCageType = "HC";
  } else if (type == cage::cageType::DoubleDiaC) {
    actualCageType = "DDC";
  } else {
    actualCageType = "error";
  }
  //
  outputFile << "1 15.999400 # " << actualCageType << "\n";
  if (padAtoms) {
    outputFile << "2 1.0 # dummy\n";
  }
  // Atoms
  outputFile << "\nAtoms\n\n";
  // -------
  // Write out the atom coordinates
  // Loop through atoms
  for (int i = 0; i < atoms.size(); i++) {
    iatom = atoms[i] - 1; // The actual index is one less than the ID
    // -----------
    // Pad out
    // Fill in dummy atoms if some have been skipped
    if (atoms[i] != prevAtomID + 1) {
      dummyAtoms = atoms[i] - prevAtomID - 1;
      dummyID = prevAtomID;
      // Loop to write out dummy atoms
      for (int j = 0; j < dummyAtoms; j++) {
        dummyID++;
        jatom = dummyID - 1;
        // 1 molecule-tag atom-type q x y z
        outputFile << dummyID << " " << yCloud->pts[jatom].molID << " 2 0 "
                   << yCloud->pts[jatom].x << " " << yCloud->pts[jatom].y << " "
                   << yCloud->pts[jatom].z << "\n";
      } // end of dummy atom write-out
    }   // end of check for dummy atom printing
    // -----------
    // Write out coordinates
    // 1 molecule-tag atom-type q x y z
    outputFile << atoms[i] << " " << yCloud->pts[iatom].molID << " 1 0 "
               << yCloud->pts[iatom].x << " " << yCloud->pts[iatom].y << " "
               << yCloud->pts[iatom].z << "\n";
    // update the previous atom ID
    prevAtomID = atoms[i];
  } // end of loop through all atoms in atomID

  // Fill in the rest of the dummy atoms
  if (atoms[atoms.size() - 1] != yCloud->nop) {
    //
    for (int id = atoms[atoms.size() - 1] + 1; id <= yCloud->nop; id++) {
      jatom = id - 1;
      outputFile << id << " " << yCloud->pts[jatom].molID << " 2 0 "
                 << yCloud->pts[jatom].x << " " << yCloud->pts[jatom].y << " "
                 << yCloud->pts[jatom].z << "\n";
    } // end of printing out dummy atoms
  }

  // Print the bonds now!
  outputFile << "\nBonds\n\n";
  // Loop through all bonds
  for (int ibond = 0; ibond < bonds.size(); ibond++) {
    // write out the bond
    outputFile << ibond + 1 << " 1 " << bonds[ibond][0] << " "
               << bonds[ibond][1] << "\n";

  } // end of for loop for bonds

  outputFile.close();
  // Once the datafile has been printed, exit
  return 0;
}

// Legacy

/**
 * @details  Function for printing out info in PairCorrel struct
 */
int sout::writeDump(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                    std::string path, std::string outFile) {
  std::ofstream outputFile;
  // ----------------
  // Make the output directory if it doesn't exist
  sout::makePath(path);
  // ----------------
  // Create a new file in the output directory
  outputFile.open(path + outFile, std::ios_base::app);

  // Append stuff
  // -----------------------
  // Header

  // The format of the LAMMPS trajectory file is:
  // ITEM: TIMESTEP
  // 0
  // ITEM: NUMBER OF ATOMS
  // 4096
  // ITEM: BOX BOUNDS pp pp pp
  // -7.9599900000000001e-01 5.0164000000000001e+01
  // -7.9599900000000001e-01 5.0164000000000001e+01
  // -7.9599900000000001e-01 5.0164000000000001e+01
  // ITEM: ATOMS id type x y z
  // 1 1 0 0 0 etc
  outputFile << "ITEM: TIMESTEP\n";
  outputFile << yCloud->currentFrame << "\n";
  outputFile << "ITEM: NUMBER OF ATOMS\n";
  outputFile << yCloud->nop << "\n";
  outputFile << "ITEM: BOX BOUNDS pp pp pp\n";
  for (int k = 0; k < yCloud->boxLow.size(); k++) {
    outputFile << yCloud->boxLow[k] << " "
               << yCloud->boxLow[k] + yCloud->box[k]; // print xlo xhi etc
    // print out the tilt factors too if it is a triclinic box
    if (yCloud->box.size() == 2 * yCloud->boxLow.size()) {
      outputFile << " "
                 << yCloud->box[k + yCloud->boxLow.size()]; // this would be +2
                                                            // for a 2D box
    }
    outputFile << "\n"; // print end of line
  }                     // end of printing box lengths
  outputFile << "ITEM: ATOMS id mol type x y z\n";
  // -----------------------
  // Atom lines
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    outputFile << yCloud->pts[iatom].atomID << " " << yCloud->pts[iatom].molID;

    // Cubic ice
    if (yCloud->pts[iatom].iceType == molSys::atom_state_type::cubic) {
      outputFile << " Ic " << yCloud->pts[iatom].x << " "
                 << yCloud->pts[iatom].y << " " << yCloud->pts[iatom].z << "\n";
    } // end of cubic ice
    // Hexagonal ice
    else if (yCloud->pts[iatom].iceType == molSys::atom_state_type::hexagonal) {
      outputFile << " Ih " << yCloud->pts[iatom].x << " "
                 << yCloud->pts[iatom].y << " " << yCloud->pts[iatom].z << "\n";
    } // end hexagonal ice
    // water
    else if (yCloud->pts[iatom].iceType == molSys::atom_state_type::water) {
      outputFile << " wat " << yCloud->pts[iatom].x << " "
                 << yCloud->pts[iatom].y << " " << yCloud->pts[iatom].z << "\n";
    } // end water
    // interfacial
    else if (yCloud->pts[iatom].iceType == molSys::atom_state_type::interfacial) {
      outputFile << " intFc " << yCloud->pts[iatom].x << " "
                 << yCloud->pts[iatom].y << " " << yCloud->pts[iatom].z << "\n";
    } // end interfacial
    // clathrate
    else if (yCloud->pts[iatom].iceType == molSys::atom_state_type::clathrate) {
      outputFile << " clathrate " << yCloud->pts[iatom].x << " "
                 << yCloud->pts[iatom].y << " " << yCloud->pts[iatom].z << "\n";
    } // end clathrate
    // interClathrate
    else if (yCloud->pts[iatom].iceType == molSys::atom_state_type::interClathrate) {
      outputFile << " interClathrate " << yCloud->pts[iatom].x << " "
                 << yCloud->pts[iatom].y << " " << yCloud->pts[iatom].z << "\n";
    } // end interClathrate
    // unclassified
    else if (yCloud->pts[iatom].iceType == molSys::atom_state_type::unclassified) {
      outputFile << " unclassified " << yCloud->pts[iatom].x << " "
                 << yCloud->pts[iatom].y << " " << yCloud->pts[iatom].z << "\n";
    } // end unclassified
    // reCubic
    else if (yCloud->pts[iatom].iceType == molSys::atom_state_type::reCubic) {
      outputFile << " reIc " << yCloud->pts[iatom].x << " "
                 << yCloud->pts[iatom].y << " " << yCloud->pts[iatom].z << "\n";
    } // end reCubic
    // reHexagonal
    else {
      outputFile << " reIh " << yCloud->pts[iatom].x << " "
                 << yCloud->pts[iatom].y << " " << yCloud->pts[iatom].z << "\n";
    } // end reHex

  } // end of loop through all atoms

  // Close the file
  outputFile.close();
  return 0;
}

/**
 * @details Function for printing out values of averaged Q6, averaged Q3 and Cij
 * values
 */
int sout::writeHisto(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                     std::vector<std::vector<int>> nList,
                     std::vector<double> avgQ6) {
  std::ofstream cijFile;
  std::ofstream q3File;
  std::ofstream q6File;
  // Create a new file in the output directory
  int nNumNeighbours;
  double avgQ3;

  cijFile.open("cij.txt", std::ofstream::out | std::ofstream::app);
  q3File.open("q3.txt", std::ofstream::out | std::ofstream::app);
  q6File.open("q6.txt", std::ofstream::out | std::ofstream::app);

  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    if (yCloud->pts[iatom].type != 1) {
      continue;
    }
    // Check for slice later
    nNumNeighbours = nList[iatom].size() - 1;
    avgQ3 = 0.0;
    for (int j = 0; j < nNumNeighbours; j++) {
      cijFile << yCloud->pts[iatom].c_ij[j].c_value << "\n";
      avgQ3 += yCloud->pts[iatom].c_ij[j].c_value;
    } // Loop through neighbours
    avgQ3 /= nNumNeighbours;
    q3File << avgQ3 << "\n";
    q6File << avgQ6[iatom] << "\n";
  } // loop through all atoms

  // Close the file
  cijFile.close();
  q3File.close();
  q6File.close();

  return 0;
}

/**
 * @details Function to print out the largest ice cluster
 */
int sout::writeCluster(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::string fileName, bool isSlice, int largestIceCluster) {
  std::ofstream clusterFile;
  // Create a new file in the output directory
  clusterFile.open(fileName, std::ofstream::out | std::ofstream::app);
  clusterFile << yCloud->currentFrame << " " << largestIceCluster << "\n";
  // Close the file
  clusterFile.close();
  return 0;
}

/**
 * @details Prints out a LAMMPS data file for all the topological bulk ice for a
 * single frame, with some default options. Only Oxygen atoms are printed out.
 * Bonds are inferred from the neighbour list
 */
int sout::writeLAMMPSdataTopoBulk(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> nList, std::vector<cage::iceType> atomTypes,
    std::string path, bool bondsBetweenDummy) {
  //
  std::ofstream outputFile;
  int iatom;            // Index, not atom ID
  int currentAtomType;  // Current atom type: a value from 1 to 4
  int numAtomTypes = 6; // DDC, HC, Mixed, dummy, mixed2 and pnc
  int bondTypes = 1;
  // Bond stuff
  std::vector<std::vector<int>> bonds; // Vector of vector, with each row
                                       // containing the atom IDs of each bond
  std::string filename =
      "system-" + std::to_string(yCloud->currentFrame) + ".data";

  // ---------------
  // Get the bonds
  if (bondsBetweenDummy) {
    bonds = bond::populateBonds(nList, yCloud);
  } // create bonds between dummy atoms
  else {
    bonds = bond::populateBonds(nList, yCloud, atomTypes);
  } // only create bonds between non-dummy atoms
  //
  // ----------------
  // Make the output directory if it doesn't exist
  sout::makePath(path);
  std::string outputDirName = path + "bulkTopo";
  sout::makePath(outputDirName);
  outputDirName = path + "bulkTopo/dataFiles/";
  sout::makePath(outputDirName);
  // ----------------
  // Write output to file inside the output directory
  outputFile.open(path + "bulkTopo/dataFiles/" + filename);
  // FORMAT:
  //  Comment Line
  //  4 atoms
  //  4 bonds
  //  0 angles
  //  0 dihedrals
  //  0 impropers
  //  1 atom types
  //  1 bond types
  //  0 angle types
  //  0 dihedral types
  //  0 improper types
  //  -1.124000 52.845002  xlo xhi
  //  0.000000 54.528999  ylo yhi
  //  1.830501 53.087501  zlo zhi

  //  Masses

  //  1 15.999400 # O

  //  Atoms

  // 1 1 1 0 20.239  1.298 6.873 # O
  // 2 1 1 0 0 5.193 6.873 # O
  // 3 1 1 0 2.249 1.298 6.873 # O

  // -------
  // Write the header
  // Write comment line
  outputFile << "Written out by D-SEAMS\n";
  // Write out the number of atoms
  outputFile << yCloud->pts.size() << " "
             << "atoms"
             << "\n";
  // Number of bonds
  outputFile << bonds.size() << " bonds"
             << "\n";
  outputFile << "0 angles\n0 dihedrals\n0 impropers\n";
  // There are maxDepth-2 total types of prisms + dummy
  outputFile << numAtomTypes << " atom types\n";
  // Bond types
  outputFile
      << bondTypes
      << " bond types\n0 angle types\n0 dihedral types\n0 improper types\n";
  // Box lengths
  outputFile << yCloud->boxLow[0] << " " << yCloud->boxLow[0] + yCloud->box[0]
             << " xlo xhi\n";
  outputFile << yCloud->boxLow[1] << " " << yCloud->boxLow[1] + yCloud->box[1]
             << " ylo yhi\n";
  outputFile << yCloud->boxLow[2] << " " << yCloud->boxLow[2] + yCloud->box[2]
             << " zlo zhi\n";
  // Masses
  outputFile << "\nMasses\n\n";
  outputFile << "1 15.999400 # dummy\n";
  outputFile << "2 15.999400 # hc \n";
  outputFile << "3 15.999400 # ddc \n";
  outputFile << "4 15.999400 # mixed \n";
  outputFile << "5 15.999400 # pnc \n";
  outputFile << "6 15.999400 # pncHexaMixed \n";
  // Atoms
  outputFile << "\nAtoms\n\n";
  // -------
  // Write out the atom coordinates
  // Loop through atoms
  for (int i = 0; i < yCloud->pts.size(); i++) {
    iatom =
        yCloud->pts[i].atomID; // The actual ID can be different from the index
    //
    // Get the atom type
    // hc atom type
    if (atomTypes[i] == cage::iceType::hc) {
      currentAtomType = 2;
    } // hc
    else if (atomTypes[i] == cage::iceType::ddc) {
      currentAtomType = 3;
    } // ddc
    // mixed
    else if (atomTypes[i] == cage::iceType::mixed) {
      currentAtomType = 4;
    } // mixed
    // pnc
    else if (atomTypes[i] == cage::iceType::pnc) {
      currentAtomType = 5;
    } // mixed
    // pnc and DDCs/HCs mixed
    else if (atomTypes[i] == cage::iceType::mixed2) {
      currentAtomType = 6;
    } // mixed
    // dummy
    else {
      currentAtomType = 1;
    } // dummy
    //
    // Write out coordinates
    // atomID molecule-tag atom-type q x y z
    outputFile << iatom << " " << yCloud->pts[i].molID << " " << currentAtomType
               << " 0 " << yCloud->pts[i].x << " " << yCloud->pts[i].y << " "
               << yCloud->pts[i].z << "\n";

  } // end of loop through all atoms in pointCloud

  // Print the bonds now!
  outputFile << "\nBonds\n\n";
  // Loop through all bonds
  for (int ibond = 0; ibond < bonds.size(); ibond++) {
    //
    outputFile << ibond + 1 << " 1 " << bonds[ibond][0] << " "
               << bonds[ibond][1] << "\n";
  } // end of for loop for bonds

  // Once the datafile has been printed, exit
  return 0;
}

/**
 * @details Function for writing out the XYZ file of a particular cluster. The
 * vector atoms contains the atom indices of the atoms to be written out
 */
int sout::writeXYZcluster(
    std::string path, molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<int> atoms, int clusterID, cage::cageType type) {
  //
  std::ofstream outputFile;
  std::string filename = "cluster-" + std::to_string(clusterID) + ".xyz";
  int nAtoms = atoms.size(); // Number of atoms
  int iatom;                 // Current atom
  // ----------------
  // Make the output directory if it doesn't exist
  sout::makePath(path);
  std::string outputDirName = path + "bulkTopo";
  sout::makePath(outputDirName);
  outputDirName = path + "bulkTopo/clusterXYZ/";
  sout::makePath(outputDirName);
  outputDirName = path + "bulkTopo/clusterXYZ/frame-" +
                  std::to_string(yCloud->currentFrame);
  sout::makePath(outputDirName);
  // ----------------
  // Write output to file inside the output directory
  outputFile.open(path + "bulkTopo/clusterXYZ/frame-" +
                  std::to_string(yCloud->currentFrame) + "/" + filename);

  // Format of an XYZ file:
  //  1970
  // generated by VMD
  //  O         43.603500       16.926201       15.215700
  //  O         39.912601       14.775100       19.379200
  outputFile << nAtoms << "\n"; // Number of atoms
  outputFile
      << "Generated by d-SEAMS. 0 type=hc and 1 type =ddc\n"; // Comment line
  //
  // Write out all the atom coordinates
  for (int i = 0; i < nAtoms; i++) {
    iatom = atoms[i]; // Atom index to be printed out
    // TODO: Should print string representation
    outputFile << static_cast<int>(type) << " " << yCloud->pts[iatom].x << " "
               << yCloud->pts[iatom].y << " " << yCloud->pts[iatom].z << "\n";
  } // end of loop through all atoms

  outputFile.close();

  return 0;
}
