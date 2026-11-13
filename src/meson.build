# d-SEAMS core library (yodaLib)

ydslib_sources = files(
  'absOrientation.cpp',
  'backward.cpp',
  'bond.cpp',
  'bop.cpp',
  'bulkTUM.cpp',
  'cage_affiliation.cpp',
  'cage_enum.cpp',
  'cluster.cpp',
  'domain.cpp',
  'franzblau.cpp',
  'gpu_residency.cpp',
  'gpu_batch.cpp',
  'generic.cpp',
  'ira_sofi.cpp',
  'local_symmetry.cpp',
  'shell_kernel.cpp',
  'cage_canon.cpp',
  'chemfiles_io.cpp',
  'chill_offload.cpp',
  'chill_plus_scheme.cpp',
  'interface_order.cpp',
  'science_front.cpp',
  'interface_wc.cpp',
  'topo_fingerprint.cpp',
  'tum_offload.cpp',
  'mol_sys.cpp',
  'sphericart_ylm.cpp',
  'ylm_cartesian.cpp',
  'neighbours.cpp',
  'nucleation_pathway.cpp',
  'cutoff_csr.cpp',
  'diamond_cna.cpp',
  'order_parameter.cpp',
  'ozaki.cpp',
  'phase.cpp',
  'pntCorrespondence.cpp',
  'density.cpp',
  'rdf.cpp',
  'rdf2d.cpp',
  'ring.cpp',
  'seams_config.cpp',
  'seams_input.cpp',
  'seams_output.cpp',
  'selection.cpp',
  'shapeMatch.cpp',
  'seams_c_api.cpp',
  'water_select.cpp',
  'site.cpp',
  'brine_front.cpp',
  'steinhardt_parallel.cpp',
  'topo_bulk.cpp',
  'cage_pack.cpp',
  'voronoi_qlm.cpp',
  'structure_desc.cpp',
  'topo_one_dim.cpp',
  'topo_two_dim.cpp',
)

_lib_link_args = _link_args
# nvc++ records the output path (src/libyodaLib.so) as DT_NEEDED.
# The Catch2 tests run with workdir input/, so that relative name
# does not resolve. Name the soname the same way gcc does.
if cpp_id == 'nvidia' or cpp_id == 'nvidia_hpc' or cpp_id == 'pgi'
  _lib_link_args += ['-Wl,-soname,libyodaLib.so']
endif

# The shell hot path has no Eigen, so AVX2 here does not fault in the allocator.
_shell_link = []
if cpp_id == 'gcc' or cpp_id == 'clang'
  _shell_avx = static_library('shell_avx', 'shell_classify.cpp',
    cpp_args: _args + ['-mavx2', '-mfma'],
    include_directories: _incdirs,
    dependencies: _deps)
  _shell_link += _shell_avx
else
  ydslib_sources += files('shell_classify.cpp')
endif

_ylm_link = []
if cpp_id == 'gcc' or cpp_id == 'clang'
  _ylm_avx = static_library('ylm_avx512', 'ylm_cartesian_avx512.cpp',
    cpp_args: _args + ['-mavx512f', '-mavx512dq', '-mfma'],
    include_directories: _incdirs,
    dependencies: _deps)
  _ylm_link += _ylm_avx
  _args += ['-DSEAMS_YLM_AVX512=1']
endif

ydslib = library('yodaLib',
  ydslib_sources,
  dependencies: _deps,
  cpp_args: _args,
  link_with: _shell_link + _ylm_link,
  link_args: _lib_link_args,
  include_directories: _incdirs,
  build_rpath: _readcondb_rpath != '' ? _readcondb_rpath : meson.current_build_dir(),
  install: true)

_linkto += ydslib

yds_dep = declare_dependency(
  link_with: [ydslib],
  include_directories: _incdirs,
  compile_args: _args,
  link_args: _link_args,
  dependencies: _deps)

# Installed interface.
#
# Building the library is not the same as shipping something another
# project can build against: without headers and a pkg-config file an
# installed copy can only be linked blind, which is why consumers vendor
# this repository as a subproject instead of depending on a packaged one.
#
# The feature defines have to travel with the headers. seams_input.hpp,
# neighbours.hpp and steinhardt_device.hpp change shape under
# SEAMS_HAS_CHEMFILES, SEAMS_HAS_READCON, SEAMS_HAS_LINKCELL and
# SEAMS_HAS_OFFLOAD, so a consumer that compiles without the set this
# library was built with sees a different declaration than the one the
# library exports. yds_dep hands them to in-tree consumers through
# compile_args; extra_cflags is the same promise to an installed one.
# The warning suppressions in _args stay behind: they are this project's
# taste, not part of its interface.
install_subdir(
  'include/internal',
  install_dir: get_option('includedir') / 'seams',
  strip_directory: true)

_public_defines = []
foreach _arg : _args
  if _arg.startswith('-DSEAMS_HAS') or _arg.startswith('-D_LIBCPP_DISABLE_AVAILABILITY')
    _public_defines += _arg
  endif
endforeach

# Naming the library object here would make meson walk the dependency
# graph into Libs.private, and that graph holds absolute paths: the
# build directory of any static subproject, and whichever environment
# supplied chemfiles and highway on the build machine. Those paths do
# not survive the install, let alone reach another machine. yodaLib is
# shared and carries its own dependencies, so the consumer needs the
# link line and nothing more.
pkgconfig = import('pkgconfig')
pkgconfig.generate(
  name: 'seams-core',
  filebase: 'seams-core',
  description: 'd-SEAMS structural elucidation engine',
  version: meson.project_version(),
  subdirs: ['seams'],
  requires: ['eigen3 >= 3.4.0'],
  libraries: ['-L${libdir}', '-lyodaLib'],
  extra_cflags: _public_defines)

# Stable C ABI. Language bindings and engine plugins link this library
# and include seams.h. They do not include the C++ headers.
# The version macros are generated from project().version, so a release
# bump cannot leave a second literal behind in the sources.
_version_parts = meson.project_version().split('.')
_version_data = configuration_data()
_version_data.set('MAJOR', _version_parts[0])
_version_data.set('MINOR', _version_parts[1])
_version_data.set('PATCH', _version_parts[2])
_version_data.set('STRING', meson.project_version())
configure_file(
  input: '../include/seams_version.h.in',
  output: 'seams_version.h',
  configuration: _version_data,
  install: true,
  install_dir: get_option('includedir'))

seams_api_lib = shared_library(
  'seams',
  'seams_api.cpp',
  dependencies: yds_dep,
  include_directories: include_directories('../include'),
  cpp_args: _args + ['-fvisibility=hidden',
                     '-I' + meson.current_build_dir()],
  link_args: _lib_link_args,
  install: true)

install_headers(files('../include/seams.h'))

pkgconfig.generate(
  seams_api_lib,
  name: 'seams',
  filebase: 'seams',
  description: 'd-SEAMS stable C API',
  version: meson.project_version(),
  extra_cflags: _public_defines)

seams_api_test = executable(
  'test_seams_api',
  files('../tests/test_seams_api.c'),
  link_with: seams_api_lib,
  include_directories: include_directories('../include'),
  c_args: ['-I' + meson.current_build_dir()])

test('seams_c_api', seams_api_test, timeout: 180)
