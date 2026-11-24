# d-SEAMS core library (yodaLib)

ydslib_sources = files(
  'absOrientation.cpp',
  'bond.cpp',
  'bop.cpp',
  'bulkTUM.cpp',
  'cage_affiliation.cpp',
  'cage_enum.cpp',
  'chordless.cpp',
  'cluster.cpp',
  'domain.cpp',
  'franzblau.cpp',
  'shell2_count.cpp',
  'lc_threads.cpp',
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
  'topo_neighbour_key.cpp',
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
  'trajectory.cpp',
)

# backward-cpp's crash handler symbolizes with backtrace_symbols, which
# musl does not provide. Without it a crash is an ordinary signal.
# It goes into the executables, never libyodaLib: its static
# SignalHandling replaces the SIGSEGV handler of whatever process loads
# the library, and Julia's GC safepoints fault on purpose.
#
# backward.hpp comes from the backward-cpp wrap (header-only interface,
# no DW/BFD/libunwind backend: none of those is enabled below, matching
# the plain execinfo.h path this project has always used). This stays
# an auto feature so a build with neither a system package nor network
# access (--wrap-mode=nodownload with nothing pre-seeded) degrades to
# "not found" instead of erroring; the packaging recipes that cannot
# reach the network (conda-forge, Nix, Yggdrasil) pre-seed this wrap's
# subproject directory from its own pinned source before configuring,
# so with_backward stays enabled there too (see recipe/build-libseams.sh,
# nix/package.nix, packaging/yggdrasil/L/libseams/build_tarballs.jl).
with_backward_opt = get_option('with_backward')
backward_has_execinfo = cppc.has_header('execinfo.h')
if with_backward_opt.enabled() and not backward_has_execinfo
  error('with_backward=enabled but execinfo.h is not available')
endif
if with_backward_opt.disabled() or not backward_has_execinfo
  backward_dep = dependency('backward-cpp-disabled-by-option', required: false)
else
  backward_dep = dependency(['backward-cpp-interface', 'Backward'],
                            required: with_backward_opt.enabled(),
                            include_type: 'system',
                            fallback: ['backward-cpp', 'backward_interface_dep'],
                            default_options: [
                              'bfd=disabled', 'dw=disabled',
                              'dwarf=disabled', 'libunwind=disabled',
                            ])
endif
backward_src = backward_dep.found() ? files('backward.cpp') : []

# The wrap's own meson.build unconditionally builds and installs a
# static libbackward-cpp (its own backward_dep, which this project
# never links: see backward_interface_dep above), so when the fallback
# subproject actually ran, that lib is an installed target too. Giving
# check_no_rebuild.py's install tests a handle to it, without invoking
# the subproject when a system package satisfied backward_dep instead
# (type_name() != 'internal') or when with_backward never ran it at
# all, keeps that check accurate instead of looking the other way.
backward_cpp_installed = []
if backward_dep.found() and backward_dep.type_name() == 'internal'
  backward_cpp_subproject = subproject('backward-cpp', required: false)
  if backward_cpp_subproject.found()
    backward_cpp_installed = [backward_cpp_subproject.get_variable('backward')]
  endif
endif

# gpu_residency.cpp calls dlopen/dlsym. A shared library that leaves
# those undefined forces every consumer to add -ldl itself.
dl_dep = cppc.find_library('dl', required: false)
if dl_dep.found()
  _deps += dl_dep
endif

_lib_link_args = _link_args
# nvc++ records the output path (src/libyodaLib.so) as DT_NEEDED.
# The Catch2 tests run with workdir input/, so that relative name
# does not resolve. Name the soname the same way gcc does.
if cpp_id == 'nvidia' or cpp_id == 'nvidia_hpc' or cpp_id == 'pgi'
  _lib_link_args += ['-Wl,-soname,libyodaLib.so']
endif

# shell_classify.cpp is built for every CPU (symm::generic) and, on x86-64,
# again as symm::avx2. That copy takes AVX2 and FMA from a target pragma
# below its includes, not from -mavx2 -mfma, so the inline and template code
# of the headers stays baseline and cannot replace the generic copies
# elsewhere in the library. symm::avx2Kernels picks the copy per call.
ydslib_sources += files('shell_classify.cpp')
_shell_link = []
_ylm_link = []
if (cpp_id == 'gcc' or cpp_id == 'clang') and host_machine.cpu_family() == 'x86_64'
  _shell_avx = static_library('shell_avx', 'shell_classify.cpp',
    cpp_args: _args + ['-DSEAMS_SHELL_AVX2_VARIANT=1'],
    include_directories: _incdirs,
    dependencies: _deps)
  _shell_link += _shell_avx
  _ylm_avx = static_library('ylm_avx512', 'ylm_cartesian_avx512.cpp',
    cpp_args: _args + ['-mavx512f', '-mavx512dq', '-mfma'],
    include_directories: _incdirs,
    dependencies: _deps)
  _ylm_link += _ylm_avx
  _args += ['-DSEAMS_YLM_AVX512=1', '-DSEAMS_SHELL_AVX2=1']
endif

# libyodaLib holds the process's one copy of linkcell and exports its lc_
# symbols. A consumer that linked the archive again carried a second copy,
# and that copy's lc_set_threads won the dynamic lookup, so the width never
# reached the searches libyodaLib runs. yds_dep therefore passes the
# linkcell headers alone.
# cutoffSolid calls ylmAdd8 only on x86-64. At -O0 a compiler keeps that
# call even when the flag in front of it is constant false, and the
# arm64 link then misses the symbol. This archive is the same file
# without SEAMS_YLM_AVX512, at -O0, so the check does not need an arm64
# machine.
if get_option('with_tests')
  _ylm_probe_args = []
  foreach _ylm_arg : _args
    if _ylm_arg != '-DSEAMS_YLM_AVX512=1'
      _ylm_probe_args += _ylm_arg
    endif
  endforeach
  ylm_o0_probe = static_library(
    'ylm_o0_probe',
    'steinhardt_parallel.cpp',
    cpp_args: _ylm_probe_args,
    override_options: ['optimization=0'],
    include_directories: _incdirs,
    dependencies: _deps,
  )
endif

_lib_origin = host_machine.system() == 'darwin' ? '@loader_path' : '$ORIGIN'
_yds_build_rpath = _readcondb_rpath != '' ? _readcondb_rpath : meson.current_build_dir()
if chemfiles_libdir != ''
  _yds_build_rpath = _yds_build_rpath + ':' + chemfiles_libdir
endif

ydslib = library('yodaLib',
  ydslib_sources,
  dependencies: _deps + [linkcell_dep],
  cpp_args: _args,
  link_with: _shell_link + _ylm_link,
  link_args: _lib_link_args,
  include_directories: _incdirs,
  build_rpath: _yds_build_rpath,
  install_rpath: _lib_origin,
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
#
# meson-uninstalled/ gets those strings unchanged, and there ${libdir} is
# a lib directory under the build root that nothing is built into. The
# -L goes through a variable that the uninstalled file sets to this build
# directory. It drops subdirs too, and takes its -I from the
# include-only dependency instead.
pkgconfig = import('pkgconfig')
pkgconfig.generate(
  name: 'seams-core',
  filebase: 'seams-core',
  description: 'd-SEAMS structural elucidation engine',
  version: meson.project_version(),
  subdirs: ['seams'],
  requires: ['eigen3 >= 3.4.0'],
  libraries: [declare_dependency(include_directories: include_directories('include/internal')),
              '-L${yodalib_libdir}', '-lyodaLib'],
  extra_cflags: _public_defines,
  variables: ['yodalib_libdir=${libdir}'],
  uninstalled_variables: [
    'yodalib_libdir=${prefix}/' + import('fs').relative_to(meson.current_build_dir(),
                                                            meson.global_build_root()),
  ])

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

_seams_map = files('libseams.map')
_seams_api_link = _lib_link_args
_seams_api_link_depends = []
if host_machine.system() not in ['darwin', 'windows'] and \
   meson.get_compiler('cpp').has_link_argument(
     '-Wl,--version-script=' + (meson.current_source_dir() / 'libseams.map'))
  _seams_api_link += ['-Wl,--version-script=' + (meson.current_source_dir() / 'libseams.map')]
  _seams_api_link_depends += _seams_map
endif
_seams_api_args = _args + ['-fvisibility=hidden', '-fvisibility-inlines-hidden',
                           '-DSEAMS_BUILDING', '-I' + meson.current_build_dir()]
seams_api_lib = shared_library(
  'seams',
  'seams_api.cpp',
  dependencies: yds_dep,
  include_directories: include_directories('../include'),
  cpp_args: _seams_api_args,
  link_args: _seams_api_link,
  link_depends: _seams_api_link_depends,
  install_rpath: _lib_origin,
  install: true)

_abidiff = find_program('abidiff', required: false)
if _abidiff.found()
  subdir('abi')
endif

install_headers(files('../include/seams.h', '../include/seams.hpp'))

pkgconfig.generate(
  seams_api_lib,
  name: 'seams',
  filebase: 'seams',
  description: 'd-SEAMS stable C API and the seams.hpp C++ wrapper',
  version: meson.project_version(),
  extra_cflags: _public_defines)

_include_ups = ''
foreach _part : get_option('libdir').split('/')
  _include_ups += '../'
endforeach
_include_ups += '../../'
_seams_lib_file = 'libseams.so'
if host_machine.system() == 'darwin'
  _seams_lib_file = 'libseams.dylib'
elif host_machine.system() == 'windows'
  _seams_lib_file = 'seams.dll'
endif
_seams_cmake = configuration_data()
_seams_cmake.set('VERSION', meson.project_version())
_seams_cmake.set('SEAMS_INCLUDE_REL', _include_ups + get_option('includedir'))
_seams_cmake.set('SEAMS_LIB_FILE', _seams_lib_file)
configure_file(
  input: '../cmake/seams-config.cmake.in',
  output: 'seams-config.cmake',
  configuration: _seams_cmake,
  install: true,
  install_dir: get_option('libdir') / 'cmake' / 'seams')
configure_file(
  input: '../cmake/seams-config-version.cmake.in',
  output: 'seams-config-version.cmake',
  configuration: _seams_cmake,
  install: true,
  install_dir: get_option('libdir') / 'cmake' / 'seams')

# Test programs. A subproject build with -Dwith_tests=false skips them.
if get_option('with_tests')
  seams_api_test = executable(
    'test_seams_api',
    files('../tests/test_seams_api.c'),
    link_with: seams_api_lib,
    include_directories: include_directories('../include'),
    c_args: ['-I' + meson.current_build_dir()])

  _capi_env = environment()
  _capi_env.set('SEAMS_SAMPLE', meson.project_source_root() / 'input' / 'traj' / 'agi_water_tiny.lammpstrj')
  _capi_env.set('SEAMS_XYZ', meson.project_source_root() / 'input' / 'xyz' / 'configuration.xyz')
  # CON frames need the readcon reader. Setting the path without it makes
  # the probe fail on a build that left the reader off.
  if readcon_dep.found()
    _capi_env.set('SEAMS_CON', meson.project_source_root() / 'input' / 'eon' / 'tiny_multi_cuh2.con')
  endif
  test('seams_c_api', seams_api_test, timeout: 180, env: _capi_env)

  # A small standalone C example: export hbond edges for a real fixture
  # frame, then confirm the quenched hook gives zero deltas for an
  # identity callback.
  capi_identity_quench_test = executable(
    'capi_identity_quench',
    files('../docs/examples/capi_identity_quench.c'),
    link_with: seams_api_lib,
    include_directories: include_directories('../include'),
    c_args: ['-I' + meson.current_build_dir()])
  test('capi_identity_quench', capi_identity_quench_test,
    args: [meson.project_source_root()])

  rings_chordless_test = executable(
    'test_rings_chordless',
    files('../tests/test_rings_chordless.c'),
    link_with: seams_api_lib,
    include_directories: include_directories('../include'),
    c_args: ['-I' + meson.current_build_dir()])
  _chordless_env = environment()
  _chordless_env.set('SEAMS_TRAJ', meson.project_source_root() / 'input' / 'traj' / 'exampleTraj.lammpstrj')
  test('seams_rings_chordless', rings_chordless_test, timeout: 120, env: _chordless_env)

  docs_ice_example = executable(
    'docs_ice_example',
    files('../docs/examples/ice.c'),
    link_with: seams_api_lib,
    include_directories: include_directories('../include'),
    c_args: ['-I' + meson.current_build_dir()],
    link_args: ['-lm'])
  test('docs_ice_example', docs_ice_example, suite: 'docs', timeout: 60)

  docs_api_cpp = executable(
    'docs_api_cpp',
    files('../docs/examples/ice.cpp'),
    link_with: seams_api_lib,
    include_directories: include_directories('../include'),
    cpp_args: ['-I' + meson.current_build_dir()],
    link_args: ['-lm'])
  test('docs_api_cpp', docs_api_cpp, suite: 'docs', timeout: 60)

  docs_api_hpp = executable(
    'docs_api_hpp',
    files('../docs/examples/cpp_api.cpp'),
    link_with: seams_api_lib,
    include_directories: include_directories('../include'),
    cpp_args: ['-I' + meson.current_build_dir()])
  test(
    'docs_api_hpp',
    docs_api_hpp,
    args: [meson.project_source_root()],
    suite: 'docs',
    timeout: 120)

  docs_quickstart_c = executable(
    'docs_quickstart_c',
    files('../docs/examples/quickstart.c'),
    link_with: seams_api_lib,
    include_directories: include_directories('../include'),
    c_args: ['-I' + meson.current_build_dir()],
    link_args: ['-lm'])
  test('docs_quickstart_c', docs_quickstart_c, suite: 'docs', timeout: 60)
endif
