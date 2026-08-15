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
  'rdf2d': 'test_rdf2d.cpp',
  'selection': 'test_selection.cpp',
  'shapeMatch': 'test_shapeMatch.cpp',
  'topo_two_dim': 'test_topo_two_dim.cpp',
  'seams_input': 'test_seams_input.cpp',
  'seams_output': 'test_seams_output.cpp',
  'bulkTUM': 'test_bulkTUM.cpp',
  'pntCorrespondence': 'test_pntCorrespondence.cpp',
  'franzblau_extra': 'test_franzblau_extra.cpp',
  'ira_sofi': 'test_ira_sofi.cpp',
  'sphericart': 'test_sphericart.cpp',
  'cage_canon': 'test_cage_canon.cpp',
  'voronoi': 'test_voronoi.cpp',
  'structure_desc': 'test_structure_desc.cpp',
  'cage_affiliation': 'test_cage_affiliation.cpp',
}

if readcon_dep.found()
  test_sources += {'readcon': 'test_readcon.cpp'}
endif

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

# --enable-new-dtags is Linux-only (macOS ld does not support it)
_test_link_args = []
if host_machine.system() == 'linux'
  _test_link_args += ['-Wl,--enable-new-dtags']
endif

foreach name, src : test_sources
  test(name,
    executable('test_' + name, src,
      dependencies: [yds_dep, catch2_dep],
      cpp_args: _args,
      link_args: _test_link_args,
      include_directories: _incdirs,
      build_rpath: '$ORIGIN/../src'),
    workdir: meson.project_source_root() / 'input',
    env: _test_env,
    timeout: 120)
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
    link_args: _test_link_args,
    include_directories: _incdirs,
    build_rpath: '$ORIGIN/../src'),
  env: _test_env)

# Performance benchmark binaries (not in regular test suite)
executable('bench_algorithms', 'bench_algorithms.cpp',
  dependencies: [yds_dep, catch2_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

executable('bench_scaling', 'bench_scaling.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

_overhead_deps = [yds_dep]
if vesin_dep.found()
  _overhead_deps += [vesin_dep]
endif
executable('bench_overhead', 'bench_overhead.cpp',
  dependencies: _overhead_deps,
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

executable('bench_cages', 'bench_cages.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

executable('bench_strong', 'bench_strong.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

executable('bench_rings_incremental', 'bench_rings_incremental.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

executable('bench_trajectory_incremental', 'bench_trajectory_incremental.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

executable('bench_lammps_io', 'bench_lammps_io.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

executable('bench_lammps_walk', 'bench_lammps_walk.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

executable('report_ice_graphs', 'report_ice_graphs.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

executable('walk_seeded', 'walk_seeded.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

executable('bench_structure_desc', 'bench_structure_desc.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

executable('compare_structure_desc', 'compare_structure_desc.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

executable('gen_figure_data', 'gen_figure_data.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')
