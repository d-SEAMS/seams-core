#include <seams_input.hpp>
#include <seams_output.hpp>

/********************************************/ /**
 *  Function for printing out ring info, when there is a
 volume slice
 Uses Boost!
 ***********************************************/
int sout::writeRingsWithSlice(std::vector<std::vector<int>> rings,
                              std::vector<bool> *flag, std::string filename) {
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
    // Skip rings that do not belong to the volume slice
    if ((*flag)[iring] == false) {
      continue;
    }

    // Otherwise, write out to the file
    for (int k = 0; k < rings[iring].size(); k++) {
      outputFile << rings[iring][k] << " ";
    } // end of loop through ring elements
    outputFile << "\n";
  } // end of loop through rings

  return 0;
}

/********************************************/ /**
 *  Function for printing out ring info, when passed a vector containing IDs of
 the rings Uses Boost!
 ***********************************************/
int sout::writeHexagonals(std::vector<std::vector<int>> rings,
                          std::vector<int> *list, std::string filename) {
  std::ofstream outputFile;
  int iring;
  int dummyValue = -10;
  // ----------------
  // Return if there are no rings
  if ((*list).size() == 0) {
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
  // Write output to file inside the output directory
  outputFile.open("../output/" + filename);

  // Format:
  // 272    214    906   1361    388      1

  for (int i = 0; i < (*list).size(); i++) {
    // Get ring index to be printed
    iring = (*list)[i];
    // Skip for rings which are actually mixed
    if (iring == dummyValue) {
      continue;
    }
    // Write iring out to the file
    for (int k = 0; k < rings[iring].size(); k++) {
      outputFile << rings[iring][k] << " ";
    } // end of loop through ring elements
    outputFile << "\n";
  } // end of loop through rings

  // Once the rings have been printed, exit
  return 0;
}

/********************************************/ /**
 *  Function for printing out ring info, when passed a vector containing IDs of
 the rings, for a volume slice Uses Boost!
 ***********************************************/
int sout::writeHexagonalsWithSlice(std::vector<std::vector<int>> rings,
                                   std::vector<bool> *flag,
                                   std::vector<int> *list,
                                   std::string filename) {
  std::ofstream outputFile;
  int iring;
  int dummyValue = -10;
  // ----------------
  // Return if there are no rings
  if ((*list).size() == 0) {
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
  // Write output to file inside the output directory
  outputFile.open("../output/" + filename);

  // Format:
  // 272    214    906   1361    388      1

  for (int i = 0; i < (*list).size(); i++) {
    // Get ring index to be printed
    iring = (*list)[i];
    if (iring == dummyValue) {
      continue;
    } // Skip dummy values
    // Check if iring is inside the slice or not
    // Skip rings that do not belong to the volume slice
    if ((*flag)[iring] == false) {
      continue;
    }
    // Write iring out to the file
    for (int k = 0; k < rings[iring].size(); k++) {
      outputFile << rings[iring][k] << " ";
    } // end of loop through ring elements
    outputFile << "\n";
  } // end of loop through rings

  // Once the rings have been printed, exit
  return 0;
}

/********************************************/ /**
                                                *  Prints out an XYZ file for
                                                *DDCs or HCs. Only Oxygen atoms
                                                *are printed out
                                                ***********************************************/
int sout::writeHexXYZ(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                      std::vector<std::vector<int>> rings,
                      std::vector<int> *list, std::string type,
                      std::string filename) {
  std::ofstream outputFile;
  std::vector<int> atoms;         // Holds all atom IDs to print
  int ringSize = rings[0].size(); // Ring size of each ring in rings
  int iring;
  int iatom; // Index, not atom ID
  // ----------------
  // Return if there are no rings
  if ((*list).size() == 0) {
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
  total_size = (*list).size() * ringSize;
  // Reserve this size inside atoms
  atoms.reserve(total_size);
  // Fill up all these atom IDs
  for (int i = 0; i < (*list).size(); i++) {
    iring = (*list)[i]; // Ring ID of ring to be appended
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
  // ----------------
  // Write output to file inside the output directory
  outputFile.open("../output/" + filename);

  //  360
  // generated by VMD
  // O         0.000000        2.597000        7.772000
  // O         0.000000        2.597000       15.094000

  // Write the header
  // Write out the number of atoms
  outputFile << atoms.size() << "\n";
  // Write comment line
  outputFile << " Written out by D-SEAMS\n";

  // Write out the atom coordinates
  // Loop through atoms
  for (int i = 0; i < atoms.size(); i++) {
    iatom = atoms[i] - 1; // The actual index is one less than the ID
    // Write out coordinates
    outputFile << " " << type << "\t" << yCloud->pts[iatom].x << "\t"
               << yCloud->pts[iatom].y << "\t" << yCloud->pts[iatom].z << "\n";
  } // end of loop through all atoms

  // Once the rings have been printed, exit
  return 0;
}

/********************************************/ /**
                                                *  Prints out a LAMMPS data
                                                *file, with some default
                                                *options. Only Oxygen atoms are
                                                *printed out
                                                ***********************************************/
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

/********************************************/ /**
 *  Function for printing out ring info, when there is no
 volume slice
 Uses Boost!
 ***********************************************/
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

/********************************************/ /**
 *  Function for printing out each info, when there is no
 volume slice
 Uses Boost!
 ***********************************************/
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
    iatomIndex = (*basal1)[iring] - 1; // C++ indices are one less
    // Write the coordinates out to the file
    outputFile << yCloud->pts[iatomIndex].x << " ";
    outputFile << yCloud->pts[iatomIndex].y << " ";
    outputFile << yCloud->pts[iatomIndex].z << " ";

    outputFile << "\n";
  } // end of loop through basal1

  // For basal 2
  for (int iring = 0; iring < ringSize; iring++) {
    iatomIndex = (*basal2)[iring] - 1; // C++ indices are one less
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
      iatomIndex = (*basal1)[iring] - 1; // C++ indices are one less
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

/********************************************/ /**
                                                *  Function for writing out all
                                                *the cages of all types to a
                                                *folder in the output
                                                ***********************************************/
int sout::writeAllCages(
    std::vector<cage::Cage> *cageList, std::vector<std::vector<int>> rings,
    std::vector<std::vector<int>> nList,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud) {
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
    if (type == cage::HexC) {
      numHC++;
      sout::writeEachCage((*cageList)[icage].rings, numHC, type, rings, yCloud);
      sout::writeBasalRingsHex((*cageList)[icage].rings, numHC, nList, rings);
    } // end of write out of HCs
    // Double diamond Cages
    else if (type == cage::DoubleDiaC) {
      numDDC++;
      sout::writeEachCage((*cageList)[icage].rings, numDDC, type, rings,
                          yCloud);
    } // end of write out of DDCs
    // Mixed Cages
    else if (type == cage::Mixed) {
      numMC++;
      sout::writeEachCage((*cageList)[icage].rings, numMC, type, rings, yCloud);
    } // end of write out of MCs
    // Error
    else {
      std::cerr << "The cage is of the wrong type\n";
      continue;
    } // some error
    // ------
  } // end of loop through all cages

  return 0;
}

/********************************************/ /**
 *  Function for printing out each cage's coordinates, when there is no
 volume slice
 Uses Boost!
 ***********************************************/
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

  if (type == cage::HexC) {
    strcpy(cageChar, "../output/cages/hexCages");
    actualCageType = "hexCages";
  } else if (type == cage::DoubleDiaC) {
    strcpy(cageChar, "../output/cages/doubleDiaCages");
    actualCageType = "doubleDiaCages";
  } else if (type == cage::Mixed) {
    strcpy(cageChar, "../output/cages/mixedCages");
    actualCageType = "mixedCages";
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
  if (type == cage::HexC) {
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

/********************************************/ /**
 *  Function for printing out the basal rings only of the hexagonal cage
 described by the number cageNum
 Uses Boost!
 ***********************************************/
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

/********************************************/ /**
 *  Function for printing out the basal rings of the prism
 blocks described by the number prismNum. Uses Boost!
 ***********************************************/
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
      auto it =
          std::find(nList[l_k - 1].begin() + 1, nList[l_k - 1].end(), m_k);

      // If the element has been found, for l1
      if (it != nList[l_k - 1].end()) {
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

  int forwardJ = (*basal2)[tempJfor] - 1;
  int backwardJ = (*basal2)[tempJback] - 1;
  int currentI = (*basal1)[iatom] - 1;

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
    std::cerr << "There is some mistake with your prism basal rings.\n";
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

/********************************************/ /**
 *  Function for printing out ring info, when there is no
 volume slice
 Uses Boost!
 ***********************************************/
int sout::writePrismNum(int nPrisms, int ringSize, std::string filename) {
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
  // # Ringsize number_of_prisms
  // 272    214    906   1361    388      1

  outputFile << "# Ringsize number_of_prisms\n";

  outputFile << ringSize << "  " << nPrisms;

  outputFile.close();

  return 0;
}

/********************************************/ /**
 *  Function for printing out bond info as pairs of atom IDs, when there is no
 volume slice
 Uses Boost!
 ***********************************************/
int sout::writeBonds(std::vector<std::vector<int>> bonds,
                     std::vector<bool> *flag, std::string filename) {
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

  // The ring vector of vector looks like:
  // 272    214    906   1361    388      1
  // The bonds corresponding to this ring are:
  // 272 214
  // 214 906
  // 1361 388
  // 388 1
  // 1 272

  // Loop through all possible bonds
  for (int ibond = 0; ibond < bonds.size(); ibond++) {
    // If the bond is a duplicate bond (flagged as false), skip it
    if ((*flag)[ibond] == false) {
      continue;
    }
    // Otherwise print the bond
    outputFile << bonds[ibond][0] << " " << bonds[ibond][1] << "\n";
  } // end of loop through bonds

  return 0;
}

/********************************************/ /**
                                                *  Prints out a LAMMPS data file
                                                *for the prisms, with some
                                                *default options. Only Oxygen
                                                *atoms are printed out
                                                ***********************************************/
int sout::writeLAMMPSdataPrisms(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud,
    std::vector<std::vector<int>> rings, bool useBondFile, std::string bondFile,
    std::vector<int> listPrism, std::string filename) {
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
    bonds = bond::populateBonds(rings);
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

  // Once the datafile has been printed, exit
  return 0;
}

/********************************************/ /**
                                                *  Prints out a LAMMPS data file
                                                *for the either the DDCs or HCs,
                                                *with some default options. Only
                                                *Oxygen atoms are printed out
                                                ***********************************************/
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
  if (type == cage::HexC) {
    actualCageType = "HC";
  } else if (type == cage::DoubleDiaC) {
    actualCageType = "DDC";
  } else if (type == cage::Mixed) {
    actualCageType = "MC";
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

  // Once the datafile has been printed, exit
  return 0;
}
/// Legacy

/********************************************/ /**
 *  Function for printing out info in PairCorrel struct
 ***********************************************/
int sout::writeDump(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                    std::string outFile) {
  std::ofstream outputFile;
  // Create a new file in the output directory
  outputFile.open(outFile, std::ios_base::app);

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
      outputFile
          << " "
          << yCloud->box[k + yCloud->boxLow
                                 .size()]; // this would be +2 for a 2D box
    }
    outputFile << "\n"; // print end of line
  }                     // end of printing box lengths
  outputFile << "ITEM: ATOMS id mol type x y z\n";
  // -----------------------
  // Atom lines
  for (int iatom = 0; iatom < yCloud->nop; iatom++) {
    outputFile << yCloud->pts[iatom].atomID << " " << yCloud->pts[iatom].molID
               << " " << yCloud->pts[iatom].iceType << " "
               << yCloud->pts[iatom].x << " " << yCloud->pts[iatom].y << " "
               << yCloud->pts[iatom].z << "\n";
  } // end of loop through all atoms

  // Close the file
  outputFile.close();
  return 0;
}

/********************************************/ /**
 *  Function for printing out values of averaged Q6, averaged Q3 and
 Cij values
 ***********************************************/
int sout::writeHisto(molSys::PointCloud<molSys::Point<double>, double> *yCloud,
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
    nNumNeighbours = yCloud->pts[iatom].neighList.size();
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

/********************************************/ /**
 * Function to print out the largest ice cluster
 ***********************************************/
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
