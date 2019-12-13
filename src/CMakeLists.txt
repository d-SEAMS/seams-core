# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# Standard cmake stuff
include(FindBLAS)
include(FindLAPACK)
include(FindBoost)
include(FindGSL)
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
  ${Eigen3_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIRS}
  )
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)
# Find packages, including conan packages
find_package(
  "Lua" REQUIRED
  "GSL" REQUIRED
  "Boost" REQUIRED COMPONENTS system filesystem
  "Eigen3 3.3" REQUIRED
  "yaml-cpp"
  "fmt"
  "rang")

link_directories(${Boost_LIBRARY_DIRS})
# Link everything
target_link_libraries(
  yodaStruct
  ${Boost_LIBRARIES}
  ${Eigen3_LIBRARIES}
  ${LUA_LIBRARIES}
  ${GSL_LIBRARIES}
  stdc++fs
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
)
# if(yaml-cpp_FOUND)
#     target_link_libraries(yodaStruct yaml-cpp)
# else()
#     message(FATAL_ERROR "yaml-cpp not found")
# endif()
