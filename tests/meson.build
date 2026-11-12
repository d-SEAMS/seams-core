# d-SEAMS tests (Catch2 v3)

catch2_dep = dependency('catch2-with-main', required: true,
  fallback: ['catch2', 'catch2_with_main_dep'])

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
  'domain': 'test_domain.cpp',
  'generic': 'test_generic.cpp',
  'order_parameter': 'test_order_parameter.cpp',
  'phase': 'test_phase.cpp',
  'rdf': 'test_rdf.cpp',
  'rdf2d': 'test_rdf2d.cpp',
  'density': 'test_density.cpp',
  'selection': 'test_selection.cpp',
  'shapeMatch': 'test_shapeMatch.cpp',
  'topo_two_dim': 'test_topo_two_dim.cpp',
  'seams_input': 'test_seams_input.cpp',
  'seams_output': 'test_seams_output.cpp',
  'bulkTUM': 'test_bulkTUM.cpp',
  'pntCorrespondence': 'test_pntCorrespondence.cpp',
  'franzblau_extra': 'test_franzblau_extra.cpp',
  'ira_sofi': 'test_ira_sofi.cpp',
  'local_symmetry': 'test_local_symmetry.cpp',
  'classify_concurrent': 'test_classify_concurrent.cpp',
  'sphericart': 'test_sphericart.cpp',
  'cage_canon': 'test_cage_canon.cpp',
  'topo_fingerprint': 'test_topo_fingerprint.cpp',
  'ion_env': 'test_ion_env.cpp',
  'guest_occupancy': 'test_guest_occupancy.cpp',
  'shell_rings': 'test_shell_rings.cpp',
  'voronoi': 'test_voronoi.cpp',
  'structure_desc': 'test_structure_desc.cpp',
  'cage_affiliation': 'test_cage_affiliation.cpp',
  'cage_enum': 'test_cage_enum.cpp',
  'c_api': 'test_c_api.cpp',
  'chill_offload': 'test_chill_offload.cpp',
  'chill_plus_rule': 'test_chill_plus_rule.cpp',
  'nucleation_pathway': 'test_nucleation_pathway.cpp',
  'thread_scale': 'test_thread_scale.cpp',
  'gpu_residency': 'test_gpu_residency.cpp',
  'config': 'test_config.cpp',
  'site': 'test_site.cpp',
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
# Coverage arc counters share a cache line. Several cores bounce that
# line and the timed tests exceed their limit. One thread keeps the
# counters and the tests inside the limit.
if get_option('b_coverage')
  _test_env.set('OMP_NUM_THREADS', '1')
endif
if csm_shape_bin != ''
  _test_env.set('CSM_SHAPE', csm_shape_bin)
endif

# --enable-new-dtags is Linux-only (macOS ld does not support it)
_test_rpath = '$ORIGIN/../src'
if _readcondb_rpath != ''
  _test_rpath = _test_rpath + ':' + _readcondb_rpath
endif
_test_link_args = _link_args
if host_machine.system() == 'linux'
  _test_link_args += ['-Wl,--enable-new-dtags']
endif
# nvc++ as the linker walks static archives once. Catch2Main needs
# symbols from Catch2; --end-group pairs with --start-group on
# LDFLAGS so ld rescans both archives.
if cpp_id == 'nvidia' or cpp_id == 'nvidia_hpc' or cpp_id == 'pgi'
  _test_link_args += ['-Wl,--end-group']
endif

foreach name, src : test_sources
  # meson allows no nested ternary; the long-running suites get more time
  _timeout = 120
  if name == 'cage_enum'
    _timeout = 300
  elif name == 'bulkTUM'
    _timeout = 300
  elif name == 'domain'
    # 16000-atom rings finish in about 15s alone. A parallel run on this
    # machine reached 97547 of 97582 assertions before the 120s limit.
    _timeout = 300
  elif name == 'cage_affiliation'
    # The seeded case passes. A parallel run was killed at 120s on the
    # last assertion.
    _timeout = 300
  elif name == 'bop' and offload_enabled
    _timeout = 600
  endif
  test(name,
    executable('test_' + name, src,
      dependencies: [yds_dep, catch2_dep],
      cpp_args: _args,
      link_args: _test_link_args,
      include_directories: _incdirs,
      build_rpath: _test_rpath),
    workdir: meson.project_source_root() / 'input',
    env: _test_env,
    timeout: _timeout)
endforeach

mpfr_dep = cppc.find_library('mpfr', required: true)
gmp_dep = cppc.find_library('gmp', required: true)
test('ozaki',
  executable('test_ozaki', 'test_ozaki.cpp',
    dependencies: [yds_dep, catch2_dep, mpfr_dep, gmp_dep],
    cpp_args: _args,
    link_args: _test_link_args,
    include_directories: _incdirs,
    build_rpath: _test_rpath),
  workdir: meson.project_source_root() / 'input',
  env: _test_env,
  timeout: 180)
executable('bench_v2', 'bench_v2.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')
executable('bench_cutoff_list', 'bench_cutoff_list.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

executable('report_symmetry_cases', 'report_symmetry_cases.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: _test_rpath)

executable('bench_shell_scale', 'bench_shell_scale.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: _test_rpath)

executable('bench_naming', 'bench_naming.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: _test_rpath)

executable('bench_roc', 'bench_roc.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: _test_rpath)

executable('report_cages', 'report_cages.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: _test_rpath)

executable('bench_robust', 'bench_robust.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: _test_rpath)

executable('report_thermal', 'report_thermal.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: _test_rpath)

executable('bench_symmetry', 'bench_symmetry.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: _test_rpath)

executable('bench_ozaki', 'bench_ozaki.cpp',
  dependencies: [yds_dep, mpfr_dep, gmp_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

if mpi_found
  _mpiexec = find_program('mpiexec', required: false)
  if _mpiexec.found()
    _domain_mpi = executable('test_domain_mpi', 'test_domain_mpi.cpp',
      dependencies: [yds_dep],
      cpp_args: _args,
      link_args: _test_link_args,
      include_directories: _incdirs,
      build_rpath: _test_rpath)
    executable('rank_bench', 'rank_bench.cpp',
      dependencies: [yds_dep],
      cpp_args: _args,
      link_args: _test_link_args,
      include_directories: _incdirs,
      build_rpath: _test_rpath)
    _mpi_env = environment({'RAYON_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1'})
    _mpi_env.prepend('LD_LIBRARY_PATH', meson.project_build_root() / 'src')
    # CI asks for domain_mpi. Meson treats that as a whole-name pattern, so
    # domain_mpi_1 never matches. Repeat the name once per rank count.
    foreach _np : ['1', '2', '4']
      test('domain_mpi', _mpiexec,
        args: ['--oversubscribe', '-n', _np, _domain_mpi],
        env: _mpi_env,
        is_parallel: false,
        timeout: 300)
    endforeach
  endif
endif

# SIMD verification test (header-only, does not need libyodaLib)
_simd_deps = [catch2_dep]
_simd_args = _args
if hwy_dep.found()
  _simd_deps += [hwy_dep]
  _simd_args += ['-DSEAMS_HAS_HWY=1']
endif
if minimage_dep.found()
  _simd_deps += [minimage_dep]
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
executable('sota_f4_zeron', 'sota_f4_zeron.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

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

executable('bench_gpu_batch', 'bench_gpu_batch.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

executable('bench_pipeline', 'bench_pipeline.cpp',
  dependencies: [yds_dep],
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

executable('bench_cage_h2h', 'bench_cage_h2h.cpp',
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

_check_invariance = executable('check_invariance', 'check_invariance.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

_fuzz_cage = executable('fuzz_cage_invariance', 'fuzz_cage_invariance.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

executable('bench_churn_sweep', 'bench_churn_sweep.cpp',
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

executable('bench_stages', 'bench_stages.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

if linkcell_dep.found()
  executable('profile_knn', 'profile_knn.cpp',
    dependencies: [yds_dep],
    cpp_args: _args,
    link_args: _test_link_args,
    include_directories: _incdirs,
    build_rpath: '$ORIGIN/../src')
endif

executable('walk_compare', 'walk_compare.cpp',
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

executable('walk_graphs', 'walk_graphs.cpp',
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

# Host cageAffiliation vs the offload cage pipeline vs a spec reference
_check_device_host = executable('check_device_host', 'check_device_host.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

# Both checks exit 1 on any mismatch. levi:<q> is the PG(2, q) incidence
# graph, with 12, 54 and 375 six-rings through each vertex for q = 2, 3, 5.
_device_host_cases = {
  'mW_cutoff': ['traj/mW_cubic.lammpstrj', '1', '1', '3.5', '2', 'cutoff'],
  'mW_knn': ['traj/mW_cubic.lammpstrj', '1', '1', '3.5', '2', 'knn'],
  'example_cutoff': ['traj/exampleTraj.lammpstrj', '1', '2', '3.5', '2', 'cutoff'],
  'example_knn': ['traj/exampleTraj.lammpstrj', '1', '2', '3.5', '2', 'knn'],
  'levi2': ['levi:2', '2'],
  'levi3': ['levi:3', '2'],
  'levi5': ['levi:5', '1'],
}
foreach name, args : _device_host_cases
  test('device_host_' + name, _check_device_host,
    args: args,
    workdir: meson.project_source_root() / 'input',
    env: _test_env,
    timeout: 300)
endforeach

# The CUDA batch kernels on the host, no device needed
_cuda_kernels = custom_target('gpu_batch_kernels',
  input: meson.project_source_root() / 'src' / 'gpu_batch.cpp',
  output: 'gpu_batch_kernels.inc',
  command: [find_program('python3'), files('extract_cuda_kernels.py'),
            '@INPUT@', '@OUTPUT@'])
test('gpu_kernels_host',
  executable('gpu_kernels_host', ['gpu_kernels_host.cpp', _cuda_kernels],
    dependencies: [yds_dep],
    cpp_args: _args,
    link_args: _test_link_args,
    include_directories: _incdirs,
    build_rpath: _test_rpath),
  env: _test_env,
  timeout: 120)

_invariance_cases = {
  'mW': ['traj/mW_cubic.lammpstrj', '1', '3.5', '1', '1'],
  'example': ['traj/exampleTraj.lammpstrj', '1', '3.5', '2', '1'],
}
foreach name, args : _invariance_cases
  test('invariance_' + name, _check_invariance,
    args: args,
    workdir: meson.project_source_root() / 'input',
    env: _test_env,
    timeout: 300)
endforeach
test('invariance_missing_input', _check_invariance,
  args: ['traj/no_such_file.lammpstrj'],
  workdir: meson.project_source_root() / 'input',
  env: _test_env,
  should_fail: true)
test('invariance_empty_frame', _check_invariance,
  args: ['traj/mW_cubic.lammpstrj', '1', '3.5', '9'],
  workdir: meson.project_source_root() / 'input',
  env: _test_env,
  should_fail: true)
# 1000 translations and 1000 permutations on every deposited frame, plus
# the eight half-box shifts. The cage stage is the cost.
test('fuzz_cage_shift_permute', _fuzz_cage,
  args: [meson.project_source_root() / 'input' / 'traj', '1000'],
  env: _test_env,
  timeout: 0)
test('neighbour_record', find_program('python3'),
  args: [files('neighbour_record.py')],
  workdir: meson.project_source_root(),
  timeout: 30)
test('docs_sources', find_program('python3'),
  args: [files('docs_sources.py')],
  workdir: meson.project_source_root(),
  timeout: 30)
test('option_names', find_program('python3'),
  args: [files('option_names.py'), meson.project_build_root() / 'seams',
         meson.project_build_root()],
  workdir: meson.project_source_root(),
  timeout: 60)
test('cli_surface', find_program('python3'),
  args: [files('cli_surface.py'), meson.project_build_root() / 'seams'],
  workdir: meson.project_source_root(),
  env: {'SEAMS_HAS_CHEMFILES': chemfiles_dep.found() ? '1' : '0'},
  timeout: 600)
test('mw215_deck', find_program('python3'),
  args: [files('test_mw215_deck.py')],
  workdir: meson.project_source_root(),
  timeout: 30)
test('mw215_frame', find_program('python3'),
  args: [files('test_mw215_deck.py'), 'frame'],
  workdir: meson.project_source_root(),
  timeout: 30)
