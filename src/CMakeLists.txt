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

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# Standard cmake stuff
include(FindBLAS)
include(FindLAPACK)
include(FindBoost)
include(FindEigen)
find_package (Eigen3 3.4 REQUIRED NO_MODULE)
add_compile_options(
  # "-Wall" "-Wpedantic" "-Wextra" "-fexceptions"
  "$<$<CONFIG:Debug>:-O0;-g3;-ggdb>"
)
# Add the new cpp files here
add_executable( yodaStruct
  main.cpp
  )

# Project Libraries
include_directories(
  include/internal
  include/external
  ${LUA_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
  )
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)
find_package(
  "Lua" REQUIRED
  "Boost" REQUIRED COMPONENTS system filesystem
  "yaml-cpp"
  "fmt")

link_directories(${Boost_LIBRARY_DIRS})
# Link everything
target_link_libraries(
  yodaStruct
  ${Boost_LIBRARIES}
  Eigen3::Eigen
  ${LUA_LIBRARIES}
  fmt
  yaml-cpp
  yodaLib
  )
add_library(yodaLib SHARED
  mol_sys.cpp
  seams_input.cpp
  seams_output.cpp
  ring.cpp
  bond.cpp
  generic.cpp
  bop.cpp
  neighbours.cpp
  opt_parser.cpp
  backward.cpp
  franzblau.cpp
  topo_one_dim.cpp
  topo_two_dim.cpp
  topo_bulk.cpp
  order_parameter.cpp
  cluster.cpp
  rdf2d.cpp
  absOrientation.cpp
  pntCorrespondence.cpp
  bulkTUM.cpp
  shapeMatch.cpp
  selection.cpp
)

install(TARGETS yodaStruct yodaLib LIBRARY DESTINATION "lib"
                      ARCHIVE DESTINATION "lib"
                      RUNTIME DESTINATION "bin"
                      COMPONENT library)

# if(yaml-cpp_FOUND)
#     target_link_libraries(yodaStruct yaml-cpp)
# else()
#     message(FATAL_ERROR "yaml-cpp not found")
# endif()
