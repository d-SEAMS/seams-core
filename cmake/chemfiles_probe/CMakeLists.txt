# Copyright (c) d-SEAMS developers
# Isolated probe. The parent must not import chemfiles before it decides
# to keep the system library: FetchContent of tag 0.10.4 also creates a
# target named chemfiles.
cmake_minimum_required(VERSION 3.16)
project(seams_chemfiles_mmcif_probe LANGUAGES CXX)

if(NOT SEAMS_CHEMFILES_MIN)
  set(SEAMS_CHEMFILES_MIN 0.10.4)
endif()
# The parent passes CMAKE_PREFIX_PATH with "|" for ";". A set
# SEAMS_CHEMFILES_DIR is searched alone.
string(REPLACE "|" ";" CMAKE_PREFIX_PATH "${SEAMS_PREFIX_PATH}")
if(SEAMS_CHEMFILES_DIR)
  find_package(chemfiles ${SEAMS_CHEMFILES_MIN} REQUIRED CONFIG
    NO_DEFAULT_PATH PATHS "${SEAMS_CHEMFILES_DIR}")
else()
  find_package(chemfiles ${SEAMS_CHEMFILES_MIN} REQUIRED CONFIG)
endif()
add_executable(seams_chemfiles_mmcif_probe
  ${CMAKE_CURRENT_LIST_DIR}/../chemfiles_mmcif_probe.cpp)
target_compile_features(seams_chemfiles_mmcif_probe PRIVATE cxx_std_17)
target_link_libraries(seams_chemfiles_mmcif_probe PRIVATE chemfiles)
set_target_properties(seams_chemfiles_mmcif_probe PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  BUILD_RPATH "$<TARGET_FILE_DIR:chemfiles>")
