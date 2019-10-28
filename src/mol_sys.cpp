#include <iostream>
#include <memory>
#include <mol_sys.hpp>

// #include <boost/filesystem.hpp>

// namespace fs = boost::filesystem;

// External Libraries

/********************************************/ /**
 *  Function for reading in a lammps file. 
 *  Reads in a specified frame (frame number and not timestep value).
 *  @param[in] filename The name of the lammps trajectory file to be read in
 *  @param[in] targetFrame The frame number whose information will be read in
 *  @param[out] yCloud The outputted PointCloud
 *  @param[in] isSlice This decides whether a slice will be created or not
 *  @param[in] coordLow Contains the lower limits of the slice, if a slice is to be created
 *  @param[in] coordHigh Contains the upper limits of the slice, if a slice is to be created
 ***********************************************/
molSys::PointCloud<molSys::Point<double>, double>
molSys::readLammpsTrj(std::string filename, int targetFrame,
                      molSys::PointCloud<molSys::Point<double>, double> *yCloud,
                      bool isSlice, std::array<double, 3> coordLow,
                      std::array<double, 3> coordHigh) {
  std::unique_ptr<std::ifstream> dumpFile;
  dumpFile = std::make_unique<std::ifstream>(filename);
  std::string line;                // Current line being read in
  std::vector<std::string> tokens; // Vector containing word tokens
  std::vector<double> numbers;     // Vector containing type double numbers
  std::vector<double> tilt;        // Vector containing tilt factors
  int currentFrame = 0;            // Current frame being read in
  int nop = -1;                    // Number of atoms in targetFrame
  bool foundFrame =
      false;            // Determines whether targetFrame has been found or not
  bool readNOP = false; // Flag for reading in the number of atoms
  bool readBox = false; // Flag for reading in the box lengths
  bool readAtoms = false; // Flag for reading in the atoms
  int xIndex, yIndex, zIndex,
      typeIndex;     // Indices for x,y,z coordinates, and LAMMPS type ID
  int molIndex = 0;  // Index for molecular ID
  int atomIndex = 0; // Index for atom ID (Only used if mol ID has not been set)
  molSys::Point<double> iPoint; // Current point being read in from the file
  xIndex = yIndex = zIndex = typeIndex = -1; // Default values
  bool isTriclinic = false; // Flag for an orthogonal or triclinic box

  if (!(molSys::file_exists(filename))) {
    std::cout
        << "Fatal Error: The file does not exist or you gave the wrong path.\n";
    // Throw exception?
    return *yCloud;
  }

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
  if (dumpFile->is_open()) {
    // ----------------------------------------------------------
    // At this point we know that the dumpfile is open
    // This loop searches for targetFrame
    while (std::getline((*dumpFile), line)) {
      // Read in lines and tokenize them
      tokens = molSys::tokenizer(line);
      // Find out which timestep number
      // you are inside
      if (tokens[0].compare("ITEM:") == 0) {
        if (tokens[1].compare("TIMESTEP") == 0) {
          // Now you are in a new timestep. Update frame number
          currentFrame++;
        }
      }

      // If targetFrame has been found
      // break out of the while loop
      if (currentFrame == targetFrame) {
        foundFrame = true;
        break; // Exit the while loop
      }
    } // End of while loop searching for targetFrame
    // ----------------------------------------------------------
    // Before filling up the PointCloud, if the vectors are filled
    // empty them
    *yCloud = molSys::clearPointCloud(yCloud);

    // ----------------------------------------------------------
    // If targetFrame has been found, read in the box lengths,
    // number of atoms and then read in atom positions, type, molID
    // By default, set molID=1 if not specified
    if (foundFrame) {
      // Run this until EOF or you reach the next timestep
      while (std::getline((*dumpFile), line)) {
        // Read in lines and tokenize them into std::string words and <double> numbers
        tokens = molSys::tokenizer(line);
        numbers = molSys::tokenizerDouble(line);

        // If you've reached the timestep line then you've reached the
        // next frame. Break out of the while loop
        if (tokens[0].compare("ITEM:") == 0) {
          if (tokens[1].compare("TIMESTEP") == 0) {
            break;
          }
        }

        // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        // Read number of particles
        if (readNOP) {
          nop = std::stoi(line.data());
          readNOP = false;
          yCloud->pts.reserve(nop);
          yCloud->nop = nop;
        }
        // Read box lengths
        if (readBox) {
          // You've reached the end of box lengths
          if (tokens[0].compare("ITEM:") == 0) {
            readBox = false;
            // If the box is triclinic, get the
            // orthogonal 'bounding box'
            if (isTriclinic) {
              // Update tilt factors
              for (int k = 0; k < tilt.size(); k++) {
                yCloud->box.push_back(tilt[k]);
              }
            } // end of check for triclinic
          }
          // Or else fill up the box lengths
          else {
            yCloud->box.push_back(numbers[1] - numbers[0]); // Update box length
            yCloud->boxLow.push_back(
                numbers[0]); // Update the lower box coordinate
            // Do this for a triclinic box only
            if (numbers.size() == 3) {
              isTriclinic = true;
              tilt.push_back(numbers[2]);
            }
          }
        }
        // Read atoms into yCloud line by line
        if (readAtoms) {
          iPoint.type = numbers[typeIndex];
          iPoint.molID = numbers[molIndex];
          iPoint.atomID = numbers[atomIndex];
          iPoint.x = numbers[xIndex];
          iPoint.y = numbers[yIndex];
          iPoint.z = numbers[zIndex];
          // Check if the particle is inside the volume Slice
          // or not
          if (isSlice) { // only if a slice has been requested
            iPoint.inSlice = molSys::atomInSlice(iPoint.x, iPoint.y, iPoint.z,
                                                 coordLow, coordHigh);
          }
          yCloud->pts.push_back(iPoint);
        }
        // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

        // Tests for reading in nop, box lengths, and atoms
        if (tokens[0].compare("ITEM:") == 0) {
          if (tokens[1].compare("NUMBER") == 0) {
            readNOP = true;
          }
        }
        if (tokens[0].compare("ITEM:") == 0) {
          if (tokens[1].compare("BOX") == 0) {
            readBox = true;
          }
        }
        if (tokens[0].compare("ITEM:") == 0) {
          if (tokens[1].compare("ATOMS") == 0) {
            readAtoms = true;
            // Now find out which index is the coordinate index etc
            for (int i = 2; i < tokens.size(); i++) {
              if (tokens[i].compare("type") == 0) {
                typeIndex = i - 2;
              }
              if (tokens[i].compare("x") == 0) {
                xIndex = i - 2;
              }
              if (tokens[i].compare("y") == 0) {
                yIndex = i - 2;
              }
              if (tokens[i].compare("z") == 0) {
                zIndex = i - 2;
              }
              if (tokens[i].compare("mol") == 0) {
                molIndex = i - 2;
              }
              if (tokens[i].compare("id") == 0) {
                atomIndex = i - 2;
              }
            } // End of for loop over tokens
            if (molIndex == 0) {
              molIndex = atomIndex;
            } // Set mol ID=atomID if not given
          }
        } // End of nested if loops for checking atom

      } // End of while
    }   // End of targetFrame found
    // ----------------------------------------------------------
  } // End of if file open statement

  // Check if you filled in the frame correctly
  if (!(foundFrame)) {
    std::cout << "You entered a frame that doesn't exist.\n";
  } // Throw exception
  if (foundFrame) {
    if (yCloud->pts.size() != yCloud->nop) {
      std::cout << "Atoms didn't get filled in properly.\n";
    }
  } // Throw exception
  yCloud->currentFrame = targetFrame;

  dumpFile->close();
  return *yCloud;
}

/********************************************/ /**
 *  Function for reading in a lammps file; and saves only the Oxygen atoms.
 This is an overloaded function. The Oxygen atom ID must be specified.
 *  @param[in] filename The name of the lammps trajectory file to be read in
 *  @param[in] targetFrame The frame number whose information will be read in
 *  @param[out] yCloud The outputted PointCloud
 *  @param[in] typeO The type ID of the Oxygen atoms
 *  @param[in] isSlice This decides whether a slice will be created or not
 *  @param[in] coordLow Contains the lower limits of the slice, if a slice is to be created
 *  @param[in] coordHigh Contains the upper limits of the slice, if a slice is to be created
 ***********************************************/
molSys::PointCloud<molSys::Point<double>, double> molSys::readLammpsTrjO(
    std::string filename, int targetFrame,
    molSys::PointCloud<molSys::Point<double>, double> *yCloud, int typeO,
    bool isSlice, std::array<double, 3> coordLow,
    std::array<double, 3> coordHigh) {
  std::unique_ptr<std::ifstream> dumpFile;
  dumpFile = std::make_unique<std::ifstream>(filename);
  std::string line;                // Current line being read in
  std::vector<std::string> tokens; // Vector containing word tokens
  std::vector<double> numbers;     // Vector containing type double numbers
  std::vector<double> tilt;        // Vector containing tilt factors
  int currentFrame = 0;            // Current frame being read in
  int nop = -1;                    // Number of atoms in targetFrame
  bool foundFrame =
      false;            // Determines whether targetFrame has been found or not
  bool readNOP = false; // Flag for reading in the number of atoms
  bool readBox = false; // Flag for reading in the box lengths
  bool readAtoms = false; // Flag for reading in the atoms
  int xIndex, yIndex, zIndex,
      typeIndex;     // Indices for x,y,z coordinates, and LAMMPS type ID
  int molIndex = 0;  // Index for molecular ID
  int atomIndex = 0; // Index for atom ID (Only used if mol ID has not been set)
  molSys::Point<double> iPoint; // Current point being read in from the file
  xIndex = yIndex = zIndex = typeIndex = -1; // Default values
  bool isTriclinic = false; // Flag for an orthogonal or triclinic box
  int nOxy = 0;             // Number of oxygen atoms

  if (!(molSys::file_exists(filename))) {
    std::cout
        << "Fatal Error: The file does not exist or you gave the wrong path.\n";
    // Throw exception?
    return *yCloud;
  }

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
  if (dumpFile->is_open()) {
    // ----------------------------------------------------------
    // At this point we know that the dumpfile is open
    // This loop searches for targetFrame
    while (std::getline((*dumpFile), line)) {
      // Read in lines and tokenize them
      tokens = molSys::tokenizer(line);
      // Find out which timestep number
      // you are inside
      if (tokens[0].compare("ITEM:") == 0) {
        if (tokens[1].compare("TIMESTEP") == 0) {
          // Now you are in a new timestep. Update frame number
          currentFrame++;
        }
      }

      // If targetFrame has been found
      // break out of the while loop
      if (currentFrame == targetFrame) {
        foundFrame = true;
        break; // Exit the while loop
      }
    } // End of while loop searching for targetFrame
    // ----------------------------------------------------------
    // Before filling up the PointCloud, if the vectors are filled
    // empty them
    *yCloud = molSys::clearPointCloud(yCloud);

    // ----------------------------------------------------------
    // If targetFrame has been found, read in the box lengths,
    // number of atoms and then read in atom positions, type, molID
    // By default, set molID=1 if not specified
    if (foundFrame) {
      // Run this until EOF or you reach the next timestep
      while (std::getline((*dumpFile), line)) {
        // Read in lines and tokenize them into std::string words and <double> numbers
        tokens = molSys::tokenizer(line);
        numbers = molSys::tokenizerDouble(line);

        // If you've reached the timestep line then you've reached the
        // next frame. Break out of the while loop
        if (tokens[0].compare("ITEM:") == 0) {
          if (tokens[1].compare("TIMESTEP") == 0) {
            break;
          }
        }

        // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        // Read number of particles
        if (readNOP) {
          nop = std::stoi(line.data());
          readNOP = false;
        }
        // Read box lengths
        if (readBox) {
          // You've reached the end of box lengths
          if (tokens[0].compare("ITEM:") == 0) {
            readBox = false;
            // If the box is triclinic, get the
            // orthogonal 'bounding box'
            if (isTriclinic) {
              // Update tilt factors
              for (int k = 0; k < tilt.size(); k++) {
                yCloud->box.push_back(tilt[k]);
              }
            } // end of check for triclinic
          }
          // Or else fill up the box lengths
          else {
            yCloud->box.push_back(numbers[1] - numbers[0]); // Update box length
            yCloud->boxLow.push_back(
                numbers[0]); // Update the lower box coordinate
            // Do this for a triclinic box only
            if (numbers.size() == 3) {
              isTriclinic = true;
              tilt.push_back(numbers[2]);
            }
          }
        }
        // Read atoms into yCloud line by line
        if (readAtoms) {
          iPoint.type = numbers[typeIndex];
          iPoint.molID = numbers[molIndex];
          iPoint.atomID = numbers[atomIndex];
          iPoint.x = numbers[xIndex];
          iPoint.y = numbers[yIndex];
          iPoint.z = numbers[zIndex];
          // Check if the particle is inside the volume Slice
          // or not
          if (isSlice) { // only if a slice has been requested
            iPoint.inSlice = molSys::atomInSlice(iPoint.x, iPoint.y, iPoint.z,
                                                 coordLow, coordHigh);
          }
          // Save only oxygen atoms
          if (iPoint.type == typeO) {
            nOxy++;
            // yCloud->pts.resize(yCloud->pts.size()+1);
            yCloud->pts.push_back(iPoint);
          }
        }
        // -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

        // Tests for reading in nop, box lengths, and atoms
        if (tokens[0].compare("ITEM:") == 0) {
          if (tokens[1].compare("NUMBER") == 0) {
            readNOP = true;
          }
        }
        if (tokens[0].compare("ITEM:") == 0) {
          if (tokens[1].compare("BOX") == 0) {
            readBox = true;
          }
        }
        if (tokens[0].compare("ITEM:") == 0) {
          if (tokens[1].compare("ATOMS") == 0) {
            readAtoms = true;
            // Now find out which index is the coordinate index etc
            for (int i = 2; i < tokens.size(); i++) {
              if (tokens[i].compare("type") == 0) {
                typeIndex = i - 2;
              }
              if (tokens[i].compare("x") == 0) {
                xIndex = i - 2;
              }
              if (tokens[i].compare("y") == 0) {
                yIndex = i - 2;
              }
              if (tokens[i].compare("z") == 0) {
                zIndex = i - 2;
              }
              if (tokens[i].compare("mol") == 0) {
                molIndex = i - 2;
              }
              if (tokens[i].compare("id") == 0) {
                atomIndex = i - 2;
              }
            } // End of for loop over tokens
            if (molIndex == 0) {
              molIndex = atomIndex;
            } // Set mol ID=atomID if not given
          }
        } // End of nested if loops for checking atom

      } // End of while
    }   // End of targetFrame found
    // ----------------------------------------------------------
  } // End of if file open statement

  // Check if you filled in the frame correctly
  if (!(foundFrame)) {
    std::cout << "You entered a frame that doesn't exist.\n";
  } // Throw exception
  if (foundFrame) {
    yCloud->nop = yCloud->pts.size();
    if (yCloud->pts.size() != nOxy) {
      std::cout << "Atoms didn't get filled in properly.\n";
    }
  } // Throw exception
  yCloud->currentFrame = targetFrame;

  dumpFile->close();
  return *yCloud;
}

/********************************************/ /**
 *  Function for clearing PointCloud if it is already 
 filled. This should be called before every frame is read in.
 *  @param[out] yCloud The cleared PointCloud
 ***********************************************/
molSys::PointCloud<molSys::Point<double>, double> molSys::clearPointCloud(
    molSys::PointCloud<molSys::Point<double>, double> *yCloud) {

  //
  std::vector<molSys::Point<double>> tempPts;
  std::vector<double> tempBox;
  //
  std::vector<double> tempBox1;

  tempPts.swap(yCloud->pts);
  tempBox.swap(yCloud->box);
  tempBox1.swap(yCloud->boxLow);

  return *yCloud;
}
