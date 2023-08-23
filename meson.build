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
subproject('eigen')
eigen_dep = dependency('eigen3',
                       version: '3.4.0',
                       required: true)
_deps += [ eigen_dep ]

fmt_dep = dependency('fmt',
                     version: '9.1.0',
                     required: true)
_deps += [ fmt_dep ]

boost_dep = dependency('boost',
  modules : [
    'filesystem',
    'system',
  ],
  required : false,
)
_deps += [ boost_dep ]

libyamlcpp_dep = dependency('yaml-cpp',
                        fallback: ['yaml-cpp',
                                   'libyamlcpp_dep'])
_deps += [ libyamlcpp_dep ]

# Include Sources

_incdirs += include_directories([ 'src/include/internal',
                               'src/include/external' ])

# --------------------- Library

ydsl_sources = [
  'src/absOrientation.cpp',
  'src/backward.cpp',
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
  'src/seams_output.cpp',
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
                           version : '5.4.6',
                           required : false)
      if lua_dep.found()
        _deps += [ lua_dep ]
      endif
      if (not lua_dep.found())
        message('Trying with subproject lua, but luafilesystem might be problematic')
        subproject('lua')
        lua_dep = dependency('lua',
                             version : '5.4.4',
                             required : false)
        if lua_dep.found()
          _deps += [ lua_dep ]
        endif
      else
        message('Lua not found, will assume spack is present')
        # Assume spack
        lua_spack_dep = cppc.find_library('lua', dirs:
             meson.source_root()+'/../.spack-env/view/lib',
                                          required: true)
        _deps += [ lua_spack_dep ]
      endif
      # Build the executable
      yds = executable('yodaStruct',
                       ['src/main.cpp'],
                       link_with : _linkto,
                       dependencies: _deps,
                       cpp_args : _args,
                       include_directories : _incdirs,
                       install : true)
    endif
endif
