# Copyright (c) d-SEAMS developers
# Optional chemfiles dependency. On by default.
# find_package uses an installed copy. Otherwise FetchContent builds
# tag 0.10.4, which carries XYZ, PDB, GRO, LAMMPS, XTC, TRR, DCD,
# Amber NetCDF, CIF and mmCIF. The 0.10.4 mmCIF writer stored a
# missing residue id as "." and the reader then parsed group_PDB.
# cmake/patches/chemfiles-0.10.4-mmcif.patch writes a numeric seq id
# and treats a missing entity id as residue 0.

option(SEAMS_WITH_CHEMFILES "Read trajectories with chemfiles" ON)

if(SEAMS_WITH_CHEMFILES)
  find_package(chemfiles QUIET CONFIG)
  if(NOT chemfiles_FOUND)
    include(FetchContent)
    set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
    FetchContent_Declare(
      chemfiles
      GIT_REPOSITORY https://github.com/chemfiles/chemfiles.git
      GIT_TAG 0.10.4
      GIT_SHALLOW TRUE
      PATCH_COMMAND sh -c "patch -p1 -N -i ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-mmcif.patch || test $? -eq 1")
    set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
    set(CHFL_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(CHFL_BUILD_DOCUMENTATION OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(chemfiles)
  endif()
endif()
