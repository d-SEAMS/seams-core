# Copyright (c) d-SEAMS developers
# Isolated probe. The parent must not import chemfiles before it decides
# to keep the system library: FetchContent of tag 0.10.4 also creates a
# target named chemfiles.
cmake_minimum_required(VERSION 3.16)
project(seams_chemfiles_mmcif_probe LANGUAGES CXX)

find_package(chemfiles REQUIRED CONFIG)
add_executable(seams_chemfiles_mmcif_probe
  ${CMAKE_CURRENT_LIST_DIR}/../chemfiles_mmcif_probe.cpp)
target_compile_features(seams_chemfiles_mmcif_probe PRIVATE cxx_std_17)
target_link_libraries(seams_chemfiles_mmcif_probe PRIVATE chemfiles)
set_target_properties(seams_chemfiles_mmcif_probe PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  BUILD_RPATH "$<TARGET_FILE_DIR:chemfiles>")
