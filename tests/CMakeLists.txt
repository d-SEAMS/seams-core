add_executable(yodaStruct_test
               main.cpp
               ${PROJECT_SOURCE_DIR}/src/molecular_system.cpp
               ${PROJECT_SOURCE_DIR}/src/molecule.cpp
               ${PROJECT_SOURCE_DIR}/src/neighbors.cpp
               ${PROJECT_SOURCE_DIR}/src/parameter.cpp
               ${PROJECT_SOURCE_DIR}/src/spherical_harmonics.cpp
               bop_correl-test.cpp
               neighbors-test.cpp
               spherical_harmonics-test.cpp
)

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
target_link_libraries(yodaStruct_test
  ${Boost_LIBRARIES}
  ${LUA_LIBRARIES}
  CONAN_PKG::blaze
  CONAN_PKG::rang
  CONAN_PKG::fmt
  CONAN_PKG::yaml-cpp
  CONAN_PKG::catch2
  yodaLib
  )


# Project Libraries
include_directories(
  ${PROJECT_SOURCE_DIR}/src/include/internal
  ${PROJECT_SOURCE_DIR}/src/include/external
  )
# Run unit tests
add_test(NAME yodaStruct-unit-tests COMMAND $<TARGET_FILE:yodaStruct_test>)
