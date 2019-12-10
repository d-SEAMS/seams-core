add_executable(yodaStruct_test
               main.cpp
               franzblau-test.cpp
               topo_one_dim-test.cpp
               topo_bulk-test.cpp
               ${PROJECT_SOURCE_DIR}/src/franzblau.cpp
               ${PROJECT_SOURCE_DIR}/src/topo_one_dim.cpp
               ${PROJECT_SOURCE_DIR}/src/topo_bulk.cpp
               ${PROJECT_SOURCE_DIR}/src/ring.cpp
               ${PROJECT_SOURCE_DIR}/src/neighbours.cpp
               ${PROJECT_SOURCE_DIR}/src/mol_sys.cpp
)

# Standard cmake stuff
# include(FindBLAS)
# include(FindLAPACK)
# include(FindBoost)
include(FindGSL)

# Find packages, including conan packages
find_package(
  "Lua" REQUIRED
  "Boost" REQUIRED
  "Catch2" REQUIRED
  "GSL" REQUIRED
  "yaml-cpp"
  "fmt"
  "rang"
  )

# Link everything
target_link_libraries(yodaStruct_test
  ${Boost_LIBRARIES}
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
  )
# Run unit tests
add_test(NAME yodaStruct-unit-tests COMMAND $<TARGET_FILE:yodaStruct_test>)
