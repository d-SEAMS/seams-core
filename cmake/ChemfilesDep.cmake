# Copyright (c) d-SEAMS developers
# Optional chemfiles dependency. On by default.
# An installed copy is kept only when it can read
# cmake/mmcif_no_entity.mmcif. chemfiles 0.10.4 wrote a missing residue
# id as "." and the reader then dereferenced the absent label_entity_id
# column. A failing probe, or no installed copy, builds tag 0.10.4 and
# applies cmake/patches/chemfiles-0.10.4-mmcif.patch. That patch writes a
# numeric seq id and treats a missing entity id as residue 0.
# SEAMS_CHEMFILES_FORCE_VENDORED skips the installed copy. The Meson
# fallback sets it after its own probe has already rejected the system
# library, so this file does not probe twice and does not import the
# rejected target before FetchContent.

option(SEAMS_WITH_CHEMFILES "Read trajectories with chemfiles" ON)
option(SEAMS_CHEMFILES_FORCE_VENDORED
  "Build patched chemfiles 0.10.4 instead of a system copy" OFF)

if(SEAMS_WITH_CHEMFILES)
  set(_seams_chemfiles_system_ok FALSE)
  if(NOT SEAMS_CHEMFILES_FORCE_VENDORED)
    set(_seams_probe_src ${CMAKE_CURRENT_LIST_DIR}/chemfiles_probe)
    set(_seams_probe_bin ${CMAKE_BINARY_DIR}/seams-chemfiles-mmcif-probe)
    set(_seams_probe_mmcif ${CMAKE_CURRENT_LIST_DIR}/mmcif_no_entity.mmcif)
    execute_process(
      COMMAND ${CMAKE_COMMAND} -S ${_seams_probe_src} -B ${_seams_probe_bin}
              -DCMAKE_BUILD_TYPE=Release
              -DCMAKE_POLICY_VERSION_MINIMUM=3.5
      RESULT_VARIABLE _seams_probe_cfg
      OUTPUT_VARIABLE _seams_probe_out
      ERROR_VARIABLE _seams_probe_err)
    if(_seams_probe_cfg EQUAL 0)
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
            "system chemfiles cannot read an mmCIF atom with no residue id (exit ${_seams_probe_run}); building patched 0.10.4")
        endif()
      endif()
    endif()
  endif()

  if(_seams_chemfiles_system_ok)
    find_package(chemfiles REQUIRED CONFIG)
  else()
    include(FetchContent)
    set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
    FetchContent_Declare(
      chemfiles
      GIT_REPOSITORY https://github.com/chemfiles/chemfiles.git
      GIT_TAG 0.10.4
      GIT_SHALLOW TRUE
      PATCH_COMMAND sh -c "patch -p1 -N -i ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-mmcif.patch || test $? -eq 1")
    set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
    set(CHFL_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(CHFL_BUILD_DOCUMENTATION OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(chemfiles)
  endif()
endif()
