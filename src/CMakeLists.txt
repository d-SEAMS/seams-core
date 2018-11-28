# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

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
  )

# Project Libraries
include_directories(
  include/internal
  include/external
  )

#
# External Libraries
#
include(Conan)
conan_cmake_run(CONANFILE ../conanfile.txt  # or relative build/conanfile.txt
                BASIC_SETUP CMAKE_TARGETS
                BUILD missing)

# YAML-CPP
find_package("yaml-cpp")
target_link_libraries(yodaStruct CONAN_PKG::yaml-cpp)

# if(yaml-cpp_FOUND)
#     target_link_libraries(yodaStruct yaml-cpp)
# else()
#     message(FATAL_ERROR "yaml-cpp not found")
# endif()

# target_link_libraries(yodaStruct ${CONAN_LIBS})

# Cotire
include(Cotire)
cotire(yodaStruct)
