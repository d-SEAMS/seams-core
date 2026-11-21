# Copyright (c) d-SEAMS developers
# Optional chemfiles dependency. On by default.
# An installed copy is kept only when cmake/chemfiles_mmcif_probe.cpp
# passes. It must register the 37 VMD molfile readers of
# cmake/patches/chemfiles-0.10.4-molfiles.patch (an unpatched 0.10.x has
# three), count the frame of cmake/molfile_water.rst7, and read
# cmake/mmcif_no_entity.mmcif. chemfiles 0.10.4 wrote a missing residue
# id as "." and the reader then dereferenced the absent label_entity_id
# column. A failing probe, or no installed copy, builds tag 0.10.4 and
# applies cmake/patches/chemfiles-0.10.4-mmcif.patch. That patch writes a
# numeric seq id and treats a missing entity id as residue 0.
# cmake/patches/chemfiles-0.10.4-mol2.patch rejects a MOL2 molecule with an
# empty counts line, which 0.10.4 indexed as counts[0].
# cmake/patches/chemfiles-0.10.4-molfiles.patch builds every molfile
# plugin that reads atoms and fixes the Molfile<F> crashes those plugins hit.
# cmake/patches/chemfiles-0.10.4-external-sources.sh then repacks the
# bundled archives: molfiles.tar.gz with the QM fields, the two log
# plugins and the hostile-file fixes (molfile-sources.patch), tng.tar.gz
# and mmtf-cpp.tar.gz with their hostile-file fixes (tng-sources.patch,
# mmtf-sources.patch). cmake/patches/chemfiles-0.10.4-hostile.patch
# stops malformed files from crashing or hanging the native readers
# (mmCIF, LAMMPS Data, XDR/XTC/TRR, TPR, DCD, Amber NetCDF, CSSR, GRO,
# SMI) and the shared text and binary file layers.
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
    set(_seams_probe_rst7 ${CMAKE_CURRENT_LIST_DIR}/molfile_water.rst7)
    # The probe is a nested cmake. It must see the same prefix and the
    # same compilers as this build, or find_package(chemfiles) looks in
    # the wrong place and the failure is silent. The prefix is quoted so
    # a semicolon-separated list stays one argument.
    execute_process(
      COMMAND ${CMAKE_COMMAND} -S ${_seams_probe_src} -B ${_seams_probe_bin}
              -DCMAKE_BUILD_TYPE=Release
              -DCMAKE_POLICY_VERSION_MINIMUM=3.5
              "-DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}"
              "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
              "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
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
          COMMAND ${_seams_probe_exe} ${_seams_probe_mmcif} ${_seams_probe_rst7}
          RESULT_VARIABLE _seams_probe_run
          OUTPUT_VARIABLE _seams_probe_out
          ERROR_VARIABLE _seams_probe_err)
        string(STRIP "${_seams_probe_out}" _seams_probe_out)
        if(_seams_probe_run EQUAL 0)
          set(_seams_chemfiles_system_ok TRUE)
        elseif(_seams_probe_run EQUAL 7)
          message(STATUS
            "system chemfiles has ${_seams_probe_out}; building patched 0.10.4, which has all of them")
        elseif(_seams_probe_run EQUAL 8)
          message(STATUS
            "system chemfiles cannot count the frames of an AMBER7 restart; building patched 0.10.4")
        else()
          message(STATUS
            "system chemfiles fails the AMBER7 restart, mmCIF, or MOL2 probe (exit ${_seams_probe_run}); building patched 0.10.4")
        endif()
      else()
        message(STATUS
          "system chemfiles probe failed to compile (exit ${_seams_probe_bld}); building patched 0.10.4")
      endif()
    else()
      # No installed chemfiles makes the nested find_package fail. That is
      # not a broken probe. Any other configure failure is.
      if(NOT _seams_probe_err MATCHES "Could not find a package configuration file provided by \"chemfiles\"")
        message(STATUS
          "system chemfiles probe failed to configure (exit ${_seams_probe_cfg}); building patched 0.10.4")
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
      PATCH_COMMAND sh -c "(patch -p1 -N -i ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-mmcif.patch || test $? -eq 1) && (patch -p1 -N -i ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-mol2.patch || test $? -eq 1) && (patch -p1 -N -i ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-frame-visit.patch || test $? -eq 1) && (patch -p1 -R -s -f --dry-run -i ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-molfiles.patch > /dev/null 2>&1 || patch -p1 -N -i ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-molfiles.patch) && sh ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-external-sources.sh molfiles molfile-sources && sh ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-external-sources.sh tng tng-sources && sh ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-external-sources.sh mmtf-cpp mmtf-sources && sh ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-external-sources.sh gemmi gemmi-sources && (patch -p1 -R -s -f --dry-run -i ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-hostile.patch > /dev/null 2>&1 || patch -p1 -N -i ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-hostile.patch) && sh ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-zlib-apple.sh && sh ${CMAKE_CURRENT_LIST_DIR}/patches/chemfiles-0.10.4-fmt-libcxx.sh")
    set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
    set(CHFL_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(CHFL_BUILD_DOCUMENTATION OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(chemfiles)
    # A conda environment already holds the unpatched libchemfiles.so.0.10,
    # and the conda compiler puts that lib directory first on the RPATH.
    # Under the upstream name the loader opens that copy, which lacks
    # seams_with_frame. The patched library needs a soname of its own.
    set_target_properties(chemfiles PROPERTIES OUTPUT_NAME seams-chemfiles)
  endif()
endif()
