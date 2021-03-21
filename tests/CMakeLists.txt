#-----------------------------------------------------------------------------------
# d-SEAMS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# A copy of the GNU General Public License is available at
# http://www.gnu.org/licenses/
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
include(FindGSL)

# Find packages
find_package(
  "Lua" REQUIRED
  "Boost" REQUIRED
  "Eigen3 3.3" REQUIRED
  "Catch2" REQUIRED
  "GSL" REQUIRED
  "yaml-cpp"
  "fmt"
  )

# Link everything
target_link_libraries(yodaStruct_test
  ${Boost_LIBRARIES}
  ${Eigen3_LIBRARIES}
  ${LUA_LIBRARIES}
  ${GSL_LIBRARIES}
  fmt
  yaml-cpp
  yodaLib
  )


# Project Libraries
include_directories(
  ${PROJECT_SOURCE_DIR}/src/include/internal
  ${PROJECT_SOURCE_DIR}/src/include/external
  ${GSL_INCLUDE_DIRS}
  ${Eigen3_INCLUDE_DIRS}
  )
# Run unit tests
add_test(NAME yodaStruct-unit-tests COMMAND $<TARGET_FILE:yodaStruct_test>)
