# Copyright (c) d-SEAMS developers
cmake_minimum_required(VERSION 3.24)
project(seams_chemfiles_dep LANGUAGES C CXX)

include(FetchContent)
include(${CMAKE_CURRENT_LIST_DIR}/ChemfilesDep.cmake)

if(TARGET chemfiles)
  get_target_property(_chemfiles_imported chemfiles IMPORTED)
  if(NOT _chemfiles_imported)
    FetchContent_GetProperties(chemfiles)
    install(TARGETS chemfiles
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib
      RUNTIME DESTINATION bin)
    if(chemfiles_SOURCE_DIR)
      install(DIRECTORY ${chemfiles_SOURCE_DIR}/include/
              DESTINATION include)
    endif()
    if(chemfiles_BINARY_DIR)
      install(FILES
              ${chemfiles_BINARY_DIR}/include/chemfiles/config.h
              ${chemfiles_BINARY_DIR}/include/chemfiles/exports.h
              DESTINATION include/chemfiles
              OPTIONAL)
    endif()
  endif()
endif()
