# Copyright (c) d-SEAMS developers
# Optional chemfiles dependency. On by default.
# An installed chemfiles is found through its CMake package, under
# SEAMS_CHEMFILES_DIR first when that is set, and must be 0.10.4 or a later
# 0.10 release. Outside a cross build it is kept only when it can read
# cmake/mmcif_no_entity.mmcif. chemfiles 0.10.4 wrote a missing residue
# id as "." and the reader then dereferenced the absent label_entity_id
# column. A cross build cannot run the probe, so the installed copy is
# accepted as found and a warning says so.
#
# SEAMS_USE_VENDORED_CHEMFILES selects the patched copy:
#   AUTO  build tag 0.10.4 with cmake/patches/chemfiles-0.10.4-*.patch
#         only when no usable installed copy exists and this is not a
#         cross build. A cross build without a usable installed copy is
#         an error, since a silent fallback is the failure it exists to
#         prevent.
#   ON    always build the patched copy.
#   OFF   never build it; a missing or unusable installed copy is an error.
# The patches write a numeric mmCIF seq id, treat a missing entity id as
# residue 0, and reject a MOL2 molecule with an empty counts line, which
# 0.10.4 indexed as counts[0].
# SEAMS_CHEMFILES_FORCE_VENDORED=ON is the older spelling of
# SEAMS_USE_VENDORED_CHEMFILES=ON. The Meson fallback sets it after its own
# probe has already rejected the system library, so this file does not
# probe twice and does not import the rejected target before FetchContent.

option(SEAMS_WITH_CHEMFILES "Read trajectories with chemfiles" ON)
option(SEAMS_CHEMFILES_FORCE_VENDORED
  "Build patched chemfiles 0.10.4 instead of a system copy" OFF)
set(SEAMS_USE_VENDORED_CHEMFILES "AUTO" CACHE STRING
  "Patched chemfiles 0.10.4: AUTO, ON or OFF")
set_property(CACHE SEAMS_USE_VENDORED_CHEMFILES PROPERTY STRINGS AUTO ON OFF)
set(SEAMS_CHEMFILES_DIR "" CACHE PATH
  "Install prefix of chemfiles, searched before the default locations")

if(SEAMS_WITH_CHEMFILES)
  string(TOUPPER "${SEAMS_USE_VENDORED_CHEMFILES}" _seams_vendored)
  if(SEAMS_CHEMFILES_FORCE_VENDORED)
    set(_seams_vendored ON)
  endif()
  if(NOT _seams_vendored MATCHES "^(AUTO|ON|OFF)$")
    message(FATAL_ERROR
      "SEAMS_USE_VENDORED_CHEMFILES must be AUTO, ON or OFF, not '${SEAMS_USE_VENDORED_CHEMFILES}'")
  endif()
  set(_seams_chemfiles_min 0.10.4)

  set(_seams_chemfiles_system_ok FALSE)
  if(NOT _seams_vendored STREQUAL "ON")
    # SEAMS_CHEMFILES_DIR is searched alone when set.
    if(SEAMS_CHEMFILES_DIR)
      set(_seams_chemfiles_search NO_DEFAULT_PATH PATHS "${SEAMS_CHEMFILES_DIR}")
    else()
      set(_seams_chemfiles_search)
    endif()
    if(CMAKE_CROSSCOMPILING)
      find_package(chemfiles ${_seams_chemfiles_min} QUIET CONFIG
        ${_seams_chemfiles_search})
      if(chemfiles_FOUND)
        message(WARNING
          "chemfiles ${chemfiles_VERSION} at ${chemfiles_DIR} is used without the mmCIF probe, which cannot run when cross-compiling")
        set(_seams_chemfiles_system_ok TRUE)
      else()
        message(STATUS "no chemfiles ${_seams_chemfiles_min} or later 0.10 release found")
      endif()
    else()
      # Only the probe project imports chemfiles until the copy is accepted:
      # FetchContent of tag 0.10.4 also creates a target named chemfiles.
      set(_seams_probe_src ${CMAKE_CURRENT_LIST_DIR}/chemfiles_probe)
      set(_seams_probe_bin ${CMAKE_BINARY_DIR}/seams-chemfiles-mmcif-probe)
      set(_seams_probe_mmcif ${CMAKE_CURRENT_LIST_DIR}/mmcif_no_entity.mmcif)
      string(REPLACE ";" "|" _seams_prefix_arg "${CMAKE_PREFIX_PATH}")
      execute_process(
        COMMAND ${CMAKE_COMMAND} -S ${_seams_probe_src} -B ${_seams_probe_bin}
                -DCMAKE_BUILD_TYPE=Release
                -DCMAKE_POLICY_VERSION_MINIMUM=3.5
                -DSEAMS_CHEMFILES_MIN=${_seams_chemfiles_min}
                "-DSEAMS_CHEMFILES_DIR=${SEAMS_CHEMFILES_DIR}"
                "-DSEAMS_PREFIX_PATH=${_seams_prefix_arg}"
        RESULT_VARIABLE _seams_probe_cfg
        OUTPUT_VARIABLE _seams_probe_out
        ERROR_VARIABLE _seams_probe_err)
      if(NOT _seams_probe_cfg EQUAL 0)
        message(STATUS "no chemfiles ${_seams_chemfiles_min} or later 0.10 release found")
      else()
        execute_process(
          COMMAND ${CMAKE_COMMAND} --build ${_seams_probe_bin}
                  --target seams_chemfiles_mmcif_probe
          RESULT_VARIABLE _seams_probe_bld
          OUTPUT_VARIABLE _seams_probe_out
          ERROR_VARIABLE _seams_probe_err)
        set(_seams_probe_exe ${_seams_probe_bin}/seams_chemfiles_mmcif_probe)
        if(NOT EXISTS ${_seams_probe_exe})
          set(_seams_probe_exe ${_seams_probe_bin}/seams_chemfiles_mmcif_probe.exe)
        endif()
        if(_seams_probe_bld EQUAL 0 AND EXISTS ${_seams_probe_exe})
          execute_process(
            COMMAND ${_seams_probe_exe} ${_seams_probe_mmcif}
            RESULT_VARIABLE _seams_probe_run
            OUTPUT_VARIABLE _seams_probe_out
            ERROR_VARIABLE _seams_probe_err)
          if(_seams_probe_run EQUAL 0)
            set(_seams_chemfiles_system_ok TRUE)
          else()
            message(STATUS
              "system chemfiles fails the mmCIF or MOL2 probe (exit ${_seams_probe_run})")
          endif()
        endif()
      endif()
      if(_seams_chemfiles_system_ok)
        find_package(chemfiles ${_seams_chemfiles_min} REQUIRED CONFIG
          ${_seams_chemfiles_search})
      endif()
    endif()
  endif()

  if(NOT _seams_chemfiles_system_ok)
    if(_seams_vendored STREQUAL "OFF")
      message(FATAL_ERROR
        "no usable chemfiles ${_seams_chemfiles_min} or later 0.10 release, and SEAMS_USE_VENDORED_CHEMFILES=OFF. "
        "Set SEAMS_CHEMFILES_DIR to its install prefix.")
    elseif(CMAKE_CROSSCOMPILING AND NOT _seams_vendored STREQUAL "ON")
      message(FATAL_ERROR
        "no chemfiles ${_seams_chemfiles_min} or later 0.10 release for the target, and a cross build does not "
        "fall back to the vendored copy on its own. Provide one through CMAKE_PREFIX_PATH or "
        "SEAMS_CHEMFILES_DIR, or pass -DSEAMS_USE_VENDORED_CHEMFILES=ON.")
    endif()
  endif()

  if(NOT _seams_chemfiles_system_ok)
    include(FetchContent)
    set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
    set(_seams_chemfiles_archive
      ${CMAKE_CURRENT_LIST_DIR}/vendor/chemfiles-0.10.4.tar.gz)
    if(EXISTS ${_seams_chemfiles_archive})
      FetchContent_Declare(
        chemfiles
        URL ${_seams_chemfiles_archive}
        URL_HASH SHA256=b8232ddaae2953538274982838aa6c2df87d300f7e2f80e92c171581e06325ba
        PATCH_COMMAND sh -c "(patch -p1 -N -i ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-mmcif.patch || test $? -eq 1) && (patch -p1 -N -i ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-mol2.patch || test $? -eq 1) && (patch -p1 -N -i ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-frame-visit.patch || test $? -eq 1) && sh ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-zlib-apple.sh && sh ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-fmt-libcxx.sh")
    else()
      FetchContent_Declare(
        chemfiles
        GIT_REPOSITORY https://github.com/chemfiles/chemfiles.git
        GIT_TAG 0.10.4
        GIT_SHALLOW TRUE
        PATCH_COMMAND sh -c "(patch -p1 -N -i ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-mmcif.patch || test $? -eq 1) && (patch -p1 -N -i ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-mol2.patch || test $? -eq 1) && (patch -p1 -N -i ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-frame-visit.patch || test $? -eq 1) && sh ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-zlib-apple.sh && sh ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-fmt-libcxx.sh")
    endif()
    set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
    set(CHFL_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(CHFL_BUILD_DOCUMENTATION OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(chemfiles)
    # A conda environment already holds the unpatched libchemfiles.so.0.10,
    # and the conda compiler puts that lib directory first on the RPATH.
    # Under the upstream name the loader opens that copy, which lacks
    # seams_with_frame. The patched library needs a soname of its own.
    # chemfiles 0.10.4 declares cmake_minimum_required(VERSION 2.8.12), so
    # MACOSX_RPATH must be enabled explicitly for @rpath/ install names on macOS.
    set_target_properties(chemfiles PROPERTIES
      OUTPUT_NAME seams-chemfiles
      MACOSX_RPATH ON)
  endif()
endif()
