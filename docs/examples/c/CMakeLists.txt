cmake_minimum_required(VERSION 3.16)
project(seams_c_family_examples C)

find_package(seams REQUIRED CONFIG)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

foreach(family neighbours rings chill_plus q4q6 cages errors)
  add_executable(${family} ${family}.c)
  target_link_libraries(${family} PRIVATE seams::seams)
  if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${family} PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()
endforeach()
