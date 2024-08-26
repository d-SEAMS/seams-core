#-----------------------------------------------------------------------------------
# d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
#
# Copyright (c) 2023--present d-SEAMS core team
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the MIT License as published by
# the Open Source Initiative.
#
# A copy of the MIT License is included in the LICENSE file of this repository.
# You should have received a copy of the MIT License along with this program.
# If not, see <https://opensource.org/licenses/MIT>.
#-----------------------------------------------------------------------------------

# The project is seams-core, the executable is yodaStruct
project('seams-core', 'cpp',
  version : '1.0.1',
  default_options : ['warning_level=3', 'cpp_std=c++14'])

host_system = host_machine.system()

# Add C++ compiler options
# Note that fast-math gives incorrect results, and the third set is for debugging
_args = ['-funroll-loops'] # ['-ffast-math', '-funroll-loops'] #['-ggdb', '-Og', '-fno-inline-functions']
_deps = [] # Dependencies
_linkto = [] # All the sub-libraries
_incdirs = [] # All the includes


# Add conditionals
if host_system == 'darwin'
  _args += ['-DOSX=TRUE']
  # Workaround for
  # error: aligned allocation function of
  # type 'void *(unsigned long, enum std::align_val_t)' is only available on macOS 10.13 or newer
  add_project_arguments(['-faligned-allocation', '-ggdb', '-Og', '-fno-inline-functions'], language: 'cpp')
endif

cppc = meson.get_compiler('cpp')

# Dependencies
# libm for Unix systems
m_dep = cppc.find_library('m', required: false)
_deps += m_dep
_deps += [declare_dependency(link_args: '-lstdc++')]

# All of which can be installed with meson wrap install <blah> originally
eigen_dep = subproject('eigen', default_options: ['default_library=static']).get_variable('eigen_dep')
_deps += [ eigen_dep ]

fmt_subproj = subproject('fmt', default_options: ['default_library=static'])
_deps += fmt_subproj.get_variable('fmt_dep')

cmake = import('cmake')
optv_cppyaml = cmake.subproject_options()
optv_cppyaml.add_cmake_defines(
    {
        'CMAKE_BUILD_TYPE': 'Release',
        'YAML_BUILD_SHARED_LIBS': 'OFF',
        'CMAKE_POSITION_INDEPENDENT_CODE': true,
        'CMAKE_CXX_FLAGS': '-w',
        'CMAKE_C_FLAGS': '-w',
    },
)
cmcppyaml = cmake.subproject('yaml-cpp', options: optv_cppyaml)
_deps += cmcppyaml.dependency('yaml-cpp')


boost_dep = dependency('boost',
  modules : [
    'filesystem',
    'system',
  ],
  required : false,
  static : true,
)
_deps += [ boost_dep ]

# Include Sources

_incdirs += include_directories([ 'src/include/internal',
                               'src/include/external' ])

# --------------------- Library

ydsl_sources = [
  'src/absOrientation.cpp',
  'src/bond.cpp',
  'src/bop.cpp',
  'src/bulkTUM.cpp',
  'src/cluster.cpp',
  'src/franzblau.cpp',
  'src/generic.cpp',
  'src/mol_sys.cpp',
  'src/neighbours.cpp',
  'src/opt_parser.cpp',
  'src/order_parameter.cpp',
  'src/pntCorrespondence.cpp',
  'src/rdf2d.cpp',
  'src/ring.cpp',
  'src/seams_input.cpp',
  'src/selection.cpp',
  'src/seams_output.cpp',
  'src/selection.cpp',
  'src/shapeMatch.cpp',
  'src/topo_bulk.cpp',
  'src/topo_one_dim.cpp',
  'src/topo_two_dim.cpp',
]

ydslib = library('yodaLib',
                ydsl_sources,
                dependencies: _deps,
                cpp_args: _args,
                include_directories : _incdirs,
                install: true
               )

yds_dep = declare_dependency(
                             link_with: [ydslib, _linkto],
                             include_directories: _incdirs,
                             compile_args: _args,
                             dependencies: _deps)
                              

# -------------------- Executable

# Only when not a subproject
if (not meson.is_subproject())
    # Lua is only used in the executable
    if get_option('with_lua')
      lua_dep = dependency('lua',
                           version : '<5.4',
                           required : true)
        _deps += [ lua_dep ]
      # Build the executable
      yds = executable('yodaStruct',
                       ['src/main.cpp',
                       'src/opt_parser.cpp'],
                       link_with : _linkto,
                       dependencies: [_deps, yds_dep],
                       cpp_args : _args,
                       include_directories : _incdirs,
                       install : true)
    endif
endif
