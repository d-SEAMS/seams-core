#-----------------------------------------------------------------------------------
# d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
#
# Copyright (c) 2018--present d-SEAMS core team
#
# SPDX-License-Identifier: MIT
#-----------------------------------------------------------------------------------

project('seams-core', 'cpp',
  version: '2.0.0',
  default_options: ['warning_level=2', 'cpp_std=c++20'],
  meson_version: '>=1.3.0')

host_system = host_machine.system()

# Global accumulator variables (eOn pattern)
_args = ['-funroll-loops', '-Wno-sign-compare', '-Wno-unused-variable',
         '-Wno-unused-parameter', '-Wno-unused-but-set-variable']
_deps = []
_linkto = []
_incdirs = []

cppc = meson.get_compiler('cpp')

# GCC vs clang warning flag differences
if cppc.get_id() == 'gcc'
  _args += ['-Wno-maybe-uninitialized']
else
  _args += ['-Wno-uninitialized']
endif

# Suppress [[nodiscard]] warnings from legacy callers (to be fixed incrementally)
_args += ['-Wno-unused-result']

# Platform conditionals
if host_system == 'darwin'
  _args += ['-DOSX=TRUE']
endif

# Core dependencies
m_dep = cppc.find_library('m', required: false)
_deps += m_dep

eigen_dep = dependency('eigen3', version: '>=3.4.0', required: true)
_deps += [eigen_dep]

blas_dep = cppc.find_library('blas', required: true)
lapack_dep = cppc.find_library('lapack', required: true)
_deps += [blas_dep, lapack_dep]

# Shared-memory parallelism. libstdc++ routes std::execution::par through TBB,
# which conda-forge does not put in this environment, so a build relying on it
# would quietly run serially. OpenMP ships with the compiler and is present.
openmp_dep = dependency('openmp', required: false)
if openmp_dep.found()
  _deps += [openmp_dep]
  _args += ['-DSEAMS_HAS_OPENMP=1']
endif

# MPI is optional. Analysis splits atoms across ranks and Allgathervs
# q_lm (Lechner-Dellago needs neighbours that may live on another rank).
# This is not LAMMPS-style spatial decomposition (Plimpton 1995).
# Default disabled so the CPU wheel never picks up a stray system MPI.
mpi_dep = dependency('mpi', language: 'cpp', required: get_option('with_mpi'))
mpi_found = mpi_dep.found()
if mpi_found
  _deps += [mpi_dep]
  _args += ['-DSEAMS_HAS_MPI=1']
endif

# OpenMP target offload. The probe must actually link the combined
# construct the library ships; a bare -fopenmp-targets flag or an
# empty `#pragma omp target {}` is not enough (missing nvptx bitcode
# still outlines to the host).
offload_opt = get_option('with_openmp_offload')
offload_enabled = false
if not offload_opt.disabled()
  if not openmp_dep.found()
    if offload_opt.enabled()
      error('with_openmp_offload=enabled but OpenMP is not available')
    endif
  else
    offload_arg = '-fopenmp-targets=' + get_option('offload_triple')
    offload_probe = '''
      #include <omp.h>
      int main() {
        const int n = 64;
        double x[64] = {};
        #pragma omp target teams distribute parallel for map(tofrom: x[0:n])
        for (int i = 0; i < n; ++i) x[i] = 1.0;
        return x[0] == 1.0 ? 0 : 1;
      }
    '''
    if cppc.links(offload_probe, args: ['-fopenmp', offload_arg],
                  name: 'OpenMP target offload link')
      _args += ['-DSEAMS_HAS_OFFLOAD=1', offload_arg]
      offload_enabled = true
    elif offload_opt.enabled()
      error('OpenMP target offload was requested but the compiler ' +
            'cannot link a target region for ' + get_option('offload_triple'))
    else
      message('OpenMP target offload probe failed for ' +
              get_option('offload_triple') + '; host OpenMP only')
    endif
  endif
endif

hwy_dep = dependency('hwy', version: '>=1.3.0', required: false)
if hwy_dep.found()
  _deps += [hwy_dep]
  _args += ['-DSEAMS_HAS_HWY=1']
endif

vesin_dep = dependency('vesin', required: false, fallback: ['vesin', 'vesin_dep'])
if vesin_dep.found()
  _deps += [vesin_dep]
  _args += ['-DSEAMS_HAS_VESIN=1']
endif

chemfiles_dep = dependency('chemfiles', required: false)
if chemfiles_dep.found()
  _deps += [chemfiles_dep]
  _args += ['-DSEAMS_HAS_CHEMFILES=1']
endif

readcon_dep = dependency('readcon-core', required: false,
                         fallback: ['readcon-core', 'readcon_dep'])
if readcon_dep.found()
  _deps += [readcon_dep]
  _args += ['-DSEAMS_HAS_READCON=1']
endif

# Report which optional backends this build actually got. Each of these
# degrades silently: without vesin the neighbour search falls back to the
# brute-force path, and without a reader the corresponding formats are simply
# absent from the library, with nothing at run time to say so.
summary({
  'SIMD distance kernel (hwy)': hwy_dep.found(),
  'Shared-memory parallelism (OpenMP)': openmp_dep.found(),
  'MPI atom split (Steinhardt)': mpi_found,
  'OpenMP target offload': offload_enabled,
  'Cell-list neighbours (vesin)': vesin_dep.found(),
  'PDB/GRO/DCD readers (chemfiles)': chemfiles_dep.found(),
  '.con reader (readcon-core)': readcon_dep.found(),
}, section: 'Optional backends', bool_yn: true)

if not readcon_dep.found()
  message('readcon-core needs cbindgen >= 0.29 on PATH, which conda-forge ' +
          'does not package; install it with `cargo install cbindgen` to ' +
          'enable the .con reader.')
endif

# Include directories
_incdirs += include_directories([
  'src/include/internal',
  'src/include/external',
])

# Library
subdir('src')

# Tests
if get_option('with_tests')
  subdir('tests')
endif

# Python bindings
if get_option('with_python')
  subdir('python')
endif
