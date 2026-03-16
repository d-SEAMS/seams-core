# d-SEAMS tests (Catch2 v3)

catch2_dep = dependency('catch2-with-main', required: true)

test_sources = {
  'franzblau': 'franzblau-test.cpp',
  'topo_bulk': 'topo_bulk-test.cpp',
  'topo_one_dim': 'topo_one_dim-test.cpp',
  'absor': 'absor-test.cpp',
  'mol_sys': 'test_mol_sys.cpp',
  'neighbours': 'test_neighbours.cpp',
  'bond': 'test_bond.cpp',
  'bop': 'test_bop.cpp',
  'cluster': 'test_cluster.cpp',
  'generic': 'test_generic.cpp',
  'order_parameter': 'test_order_parameter.cpp',
}

# Ensure test binaries load the freshly-built libyodaLib.so instead of any
# version installed in the conda env.  The conda toolchain injects
# --disable-new-dtags + -rpath=<pixi-env>/lib via LDFLAGS, which creates a
# DT_RPATH entry that overrides LD_LIBRARY_PATH.  We counter this by:
#   1. Adding --enable-new-dtags so the linker emits DT_RUNPATH (which
#      LD_LIBRARY_PATH can override).
#   2. Prepending the build-dir src/ to LD_LIBRARY_PATH in the test env.
#   3. Setting build_rpath so the binary also works when run outside meson test.
_test_env = environment()
_test_env.prepend('LD_LIBRARY_PATH', meson.project_build_root() / 'src')

foreach name, src : test_sources
  test(name,
    executable('test_' + name, src,
      dependencies: [yds_dep, catch2_dep],
      cpp_args: _args,
      link_args: ['-Wl,--enable-new-dtags'],
      include_directories: _incdirs,
      build_rpath: '$ORIGIN/../src'),
    workdir: meson.project_source_root() / 'input',
    env: _test_env)
endforeach

# SIMD verification test (header-only, does not need libyodaLib)
_simd_deps = [catch2_dep]
_simd_args = _args
if hwy_dep.found()
  _simd_deps += [hwy_dep]
  _simd_args += ['-DSEAMS_HAS_HWY=1']
endif
test('simd',
  executable('test_simd', 'test_simd.cpp',
    dependencies: _simd_deps,
    cpp_args: _simd_args,
    link_args: ['-Wl,--enable-new-dtags'],
    include_directories: _incdirs,
    build_rpath: '$ORIGIN/../src'),
  env: _test_env)

# Performance benchmark binary (not in regular test suite)
executable('bench_algorithms', 'bench_algorithms.cpp',
  dependencies: [yds_dep, catch2_dep],
  cpp_args: _args,
  link_args: ['-Wl,--enable-new-dtags'],
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')
