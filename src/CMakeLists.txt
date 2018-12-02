# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# Standard cmake stuff
include(FindBLAS)
include(FindLAPACK)
include(FindBoost)

# Add the new cpp files here
add_executable( yodaStruct
  main.cpp
  spherical_harmonics.cpp
  )

# Project Libraries
include_directories(
  include/internal
  include/external
  ${LUA_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
  )

# Conan Stuff
include(Conan)
conan_cmake_run(CONANFILE ../conanfile.txt  # or relative build/conanfile.txt
                BASIC_SETUP CMAKE_TARGETS
                BUILD missing)

# Find packages, including conan packages
find_package(
  "Lua" REQUIRED
  "Boost" REQUIRED
  "blaze"
  "yaml-cpp"
  "fmt"
  "rang")

# Link everything
target_link_libraries(
  yodaStruct
  ${Boost_LIBRARIES}
  ${LUA_LIBRARIES}
  CONAN_PKG::blaze
  CONAN_PKG::rang
  CONAN_PKG::fmt
  CONAN_PKG::yaml-cpp
  yodaLib
  )

add_library(yodaLib SHARED
  bop_correl.cpp
  molecular_system.cpp
  molecule.cpp
  neighbors.cpp
  parameter.cpp
  structure_factor.cpp
  transition.cpp
  opt_parser.cpp
  output.cpp
  rdf2D.cpp
  rdf3D.cpp
  backward.cpp
  density.cpp
  geometry.cpp
)
# if(yaml-cpp_FOUND)
#     target_link_libraries(yodaStruct yaml-cpp)
# else()
#     message(FATAL_ERROR "yaml-cpp not found")
# endif()

# target_link_libraries(yodaStruct ${CONAN_LIBS})
