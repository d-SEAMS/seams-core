project('yodaStruct', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++17'])

host_system = host_machine.system()
# Note that fast-math gives incorrect results, and the third set is for debugging
yoda_extra_args = ['-funroll-loops'] # ['-ffast-math', '-funroll-loops'] #['-ggdb', '-Og', '-fno-inline-functions']

add_languages('c', required: true)

cc = meson.get_compiler('c')
cppc = meson.get_compiler('cpp')

# Add conditionals
if host_system == 'darwin'
  yoda_extra_args += ['-DOSX=TRUE']
  # Workaround for
  # error: aligned allocation function of
  # type 'void *(unsigned long, enum std::align_val_t)' is only available on macOS 10.13 or newer
  add_global_arguments(['-faligned-allocation', '-ggdb', '-Og', '-fno-inline-functions'], language: 'cpp')
endif

# Dependencies
eigen_dep = dependency('Eigen3', modules : ['Eigen3::Eigen'])
yds_deps = [ eigen_dep ]

boost_dep = dependency('boost',
  modules : [
    'filesystem',
    'system',
  ],
  required : false,
)

fmt_dep = dependency('fmt', method : 'cmake',  modules : ['fmt::fmt'],  required : true)

libyamlcpp = dependency('yaml-cpp',
                        fallback: ['yaml-cpp', 'libyamlcpp_dep'])

lua_dep = dependency('lua',
                     version : '<5.3',
                     required : false)
if lua_dep.found()
  yds_deps += [ lua_dep ]
else
  # Assume spack
  lua_spack_dep = cppc.find_library('lua', dirs:
       meson.source_root()+'/../.spack-env/view/lib',
                                    required: true)
  yds_deps += [lua_spack_dep]
endif

yds_deps += [ boost_dep, fmt_dep, libyamlcpp ]

incdir = include_directories([ 'include/internal', 'include/external' ])

ydsl_sources = [
'absOrientation.cpp',
'backward.cpp',
'bond.cpp',
'bop.cpp',
'bulkTUM.cpp',
'cluster.cpp',
'franzblau.cpp',
'generic.cpp',
'mol_sys.cpp',
'neighbours.cpp',
'opt_parser.cpp',
'order_parameter.cpp',
'pntCorrespondence.cpp',
'rdf2d.cpp',
'ring.cpp',
'seams_input.cpp',
'seams_output.cpp',
'shapeMatch.cpp',
'topo_bulk.cpp',
'topo_one_dim.cpp',
'topo_two_dim.cpp',
'bulkClathrate.cpp',
'selection.cpp'
]

# ---------------------- Executable

ydslib = library('yodaLib',
                ydsl_sources,
                dependencies: yds_deps,
                cpp_args: yoda_extra_args,
                include_directories : incdir,
                install: true
               )

yds = executable('yodaStruct',
                ['main.cpp'],
                link_with : ydslib,
                dependencies: yds_deps,
                cpp_args : yoda_extra_args,
                include_directories : incdir,
                install : true)

