/// C++ Structure Determination Toolkit
///
/// Copyright © 2018-2019 Rohit Goswami <r95g10[at]gmail.com>, Amrita Goswami
/// <amritag[at]iitk.ac.in>
///
/// Main Catch test driver
///
/// @file main.cpp
/// @brief Catch test driver
/// @author Rohit Goswami
#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#define CATCH_CONFIG_RUNNER
#define CATCH_CONFIG_NO_CPP17_UNCAUGHT_EXCEPTIONS
#define CATCH_CONFIG_CPP17_STRING_VIEW
#include <catch2/catch.hpp>

// int main(int argc, char* argv[]) {
//   Catch::Session session;  // There must be exactly one instance

//   // writing to session.configData() here sets defaults
//   // this is the preferred way to set them

//   int returnCode = session.applyCommandLine(argc, argv);
//   if (returnCode != 0) {  // Indicates a command line error
//     return returnCode;
//   }

//   // writing to session.configData() or session.Config() here
//   // overrides command line args
//   // only do this if you know you need to

//   int numFailed = session.run();

//   // numFailed is clamped to 255 as some unixes only use the lower 8 bits.
//   // This clamping has already been applied, so just return it here
//   // You can also do any post run clean-up here
//   return numFailed;
// }
