# Find Eigen library
# This module sets the following variables:
#   EIGEN_FOUND        : True if Eigen library is found, False otherwise
#   EIGEN_INCLUDE_DIRS : Include directories for Eigen headers
#   EIGEN_LIBRARIES    : Eigen library to link against

# Search for Eigen headers
find_path(EIGEN_INCLUDE_DIR NAMES Eigen/Core PATH_SUFFIXES eigen3)

# Check if Eigen headers are found
if(EIGEN_INCLUDE_DIR)
  set(EIGEN_FOUND TRUE)
  set(EIGEN_INCLUDE_DIRS ${EIGEN_INCLUDE_DIR})
else()
  set(EIGEN_FOUND FALSE)
  set(EIGEN_INCLUDE_DIRS)
endif()

# Provide the information to CMake
if(EIGEN_FOUND)
  message(STATUS "Found Eigen: ${EIGEN_INCLUDE_DIRS}")
else()
  message(STATUS "Eigen not found. Please set EIGEN_INCLUDE_DIRS manually.")
endif()

# Export the variables
set(EIGEN_INCLUDE_DIRS ${EIGEN_INCLUDE_DIRS} CACHE PATH "Eigen include directories" FORCE)
mark_as_advanced(EIGEN_INCLUDE_DIRS)

