# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# Add the new cpp files here
add_executable( yodaStruct
  backward.cpp
  main.cpp
  molecule.cpp
  molecular_system.cpp
  parameter.cpp
  rdf3D.cpp
  rdf2D.cpp
  structure_factor.cpp
  opt_parser.cpp
  output.cpp
  geometry.cpp
  density.cpp
  transition.cpp
  )

# Project Libraries
include_directories(
  include/internal
  include/external
  ${LUA_INCLUDE_DIR}
  )

# Conan Stuff
include(Conan)
conan_cmake_run(CONANFILE ../conanfile.txt  # or relative build/conanfile.txt
                BASIC_SETUP CMAKE_TARGETS
                BUILD missing)

# Find packages, including conan packages
find_package(
  "Lua" REQUIRED
  "yaml-cpp"
  "fmt"
  "rang")

# Link everything
target_link_libraries(yodaStruct
  ${LUA_LIBRARIES}
  CONAN_PKG::rang
  CONAN_PKG::fmt
  CONAN_PKG::yaml-cpp)

# if(yaml-cpp_FOUND)
#     target_link_libraries(yodaStruct yaml-cpp)
# else()
#     message(FATAL_ERROR "yaml-cpp not found")
# endif()

# target_link_libraries(yodaStruct ${CONAN_LIBS})
