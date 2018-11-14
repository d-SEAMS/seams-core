# Add the new cpp files here
add_executable( yodaStruct
  main.cpp
  molecule.cpp
  molecular_system.cpp
  parameter.cpp
  rdf3D.cpp
  rdf2D.cpp
  structure_factor.cpp
  output.cpp
  geometry.cpp
  density.cpp
  transition.cpp
  chunk.cpp
  )

# # Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

include(Conan)
# conan_cmake_run(REQUIRES Hello/0.1@memsharded/testing
#                 BASIC_SETUP
#                 BUILD missing)
