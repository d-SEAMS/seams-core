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

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
