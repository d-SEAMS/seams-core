#-----------------------------------------------------------------------------------
# d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
#
# Copyright (c) 2018--present d-SEAMS core team
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the MIT License as published by
# the Open Source Initiative.
#
# A copy of the MIT License is included in the LICENSE file of this repository.
# You should have received a copy of the MIT License along with this program.
# If not, see <https://opensource.org/licenses/MIT>.
#-----------------------------------------------------------------------------------

add_executable(yodaStruct_test
               main.cpp
               franzblau-test.cpp
               topo_one_dim-test.cpp
               topo_bulk-test.cpp
               absor-test.cpp
               ${PROJECT_SOURCE_DIR}/src/franzblau.cpp
               ${PROJECT_SOURCE_DIR}/src/topo_one_dim.cpp
               ${PROJECT_SOURCE_DIR}/src/topo_bulk.cpp
               ${PROJECT_SOURCE_DIR}/src/ring.cpp
               ${PROJECT_SOURCE_DIR}/src/neighbours.cpp
               ${PROJECT_SOURCE_DIR}/src/mol_sys.cpp
               ${PROJECT_SOURCE_DIR}/src/absOrientation.cpp
               ${PROJECT_SOURCE_DIR}/src/seams_input.cpp
               ${PROJECT_SOURCE_DIR}/src/seams_output.cpp
               ${PROJECT_SOURCE_DIR}/src/pntCorrespondence.cpp
               ${PROJECT_SOURCE_DIR}/src/bulkTUM.cpp
)

# Standard cmake stuff
# include(FindBLAS)
# include(FindLAPACK)
# include(FindBoost)

# Find packages
find_package(
  "Lua" REQUIRED
  "Boost" REQUIRED
  "Eigen3 3.3" REQUIRED
  "Catch2" REQUIRED
  "yaml-cpp"
  "fmt"
  )

# Link everything
target_link_libraries(yodaStruct_test
  ${Boost_LIBRARIES}
  ${Eigen3_LIBRARIES}
  ${LUA_LIBRARIES}
  fmt
  yaml-cpp
  yodaLib
  )


# Project Libraries
include_directories(
  ${PROJECT_SOURCE_DIR}/src/include/internal
  ${PROJECT_SOURCE_DIR}/src/include/external
  ${Eigen3_INCLUDE_DIRS}
  )
# Run unit tests
add_test(NAME yodaStruct-unit-tests COMMAND $<TARGET_FILE:yodaStruct_test>)
