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

readcon_dep = dependency('readcon-core', required: false)
if readcon_dep.found()
  _deps += [readcon_dep]
  _args += ['-DSEAMS_HAS_READCON=1']
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
