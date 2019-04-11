add_executable(yodaStruct_test
               main.cpp
               bop_correl-test.cpp
               neighbors-test.cpp
               spherical_harmonics-test.cpp
               ${PROJECT_SOURCE_DIR}/src/spherical_harmonics.cpp
)

# Find packages
find_package(
  "Lua" REQUIRED
  "Boost" REQUIRED
  "yaml-cpp"
  "fmt"
  "rang")

# Link everything
target_link_libraries(yodaStruct_test
  ${Boost_LIBRARIES}
  ${LUA_LIBRARIES}
  fmt
  yaml-cpp
  catch2
  rang
  yodaLib
  )


# Project Libraries
include_directories(
  ${PROJECT_SOURCE_DIR}/src/include/internal
  ${PROJECT_SOURCE_DIR}/src/include/external
  )
# Run unit tests
add_test(NAME yodaStruct-unit-tests COMMAND $<TARGET_FILE:yodaStruct_test>)
