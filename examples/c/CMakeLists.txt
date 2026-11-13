cmake_minimum_required(VERSION 3.16)
project(seams_c_quickstart C)

find_package(seams REQUIRED)

add_executable(quickstart quickstart.c)
target_link_libraries(quickstart PRIVATE seams::seams m)
