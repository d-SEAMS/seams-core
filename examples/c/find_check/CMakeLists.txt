cmake_minimum_required(VERSION 3.16)
project(seams_find_check C)

if(NOT DEFINED SEAMS_EXPECT_VERSION)
  message(FATAL_ERROR "pass -DSEAMS_EXPECT_VERSION=<major.minor.patch>")
endif()
find_package(seams ${SEAMS_EXPECT_VERSION} EXACT REQUIRED)
if(NOT TARGET seams::seams)
  message(FATAL_ERROR "seams-config.cmake did not define seams::seams")
endif()

add_executable(version version.c)
target_link_libraries(version PRIVATE seams::seams)
