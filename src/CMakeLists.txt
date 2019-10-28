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
  )

# Project Libraries
include_directories(
  include/internal
  include/external
  ${LUA_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
  )

# Find packages, including conan packages
find_package(
  "Lua" REQUIRED
  "Boost" REQUIRED
  "yaml-cpp"
  "fmt"
  "rang")

# Link everything
target_link_libraries(
  yodaStruct
  ${Boost_LIBRARIES}
  ${LUA_LIBRARIES}
  fmt
  yaml-cpp
  yodaLib
  )
add_library(yodaLib SHARED
  mol_sys.cpp
  seams-input.cpp
  generic.cpp
  bop.cpp
  neighbours.cpp
  opt_parser.cpp
  backward.cpp
)
# if(yaml-cpp_FOUND)
#     target_link_libraries(yodaStruct yaml-cpp)
# else()
#     message(FATAL_ERROR "yaml-cpp not found")
# endif()
