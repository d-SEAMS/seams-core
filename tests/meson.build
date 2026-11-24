# d-SEAMS tests (Catch2 v3)

# Prefer an installed Catch2 that actually links. A copy built against a
# newer libc, or without PIE, is replaced by the wrap, which only the tests
# use and so is not installed with d-SEAMS.
catch2_dep = dependency('catch2-with-main', required: false)
if not catch2_dep.found() or not cppc.links(
    '#include <catch2/catch_test_macros.hpp>\nTEST_CASE("link") {}\n',
    dependencies: catch2_dep,
    name: 'installed catch2 links')
  catch2_dep = subproject('catch2', default_options: ['tests=false', 'install=false']).get_variable(
    'catch2_with_main_dep')
endif

# rang.hpp colourizes franzblau/absor/topo_bulk/topo_one_dim's test output.
rang_dep = dependency('rang', required: true, include_type: 'system',
                      fallback: ['rang', 'rang_dep'])
_rang_names = ['franzblau', 'topo_bulk', 'topo_one_dim', 'absor']

# cxxopts.hpp backs tests/include/internal/options.h's Argv helper. Wrapped
# so a future test that includes it needs no plumbing beyond this list; the
# CLI itself parses with Argum (third_party/argum/argum.h), not cxxopts.
cxxopts_dep = dependency('cxxopts', required: true, include_type: 'system',
                         fallback: ['cxxopts', 'cxxopts_dep'])
# Compiled at configure time so --wrap-mode=forcefallback and a plain
# system cxxopts are both exercised, even though no test TU includes
# tests/include/internal/options.h right now.
if not cppc.compiles(
    '#include <cxxopts.hpp>\n' +
    'int main(int argc, char **argv) {\n' +
    '  cxxopts::Options opts("seams-tests", "cxxopts wrap smoke check");\n' +
    '  opts.add_options()("h,help", "Print usage");\n' +
    '  return opts.parse(argc, argv).count("help") >= 0 ? 0 : 1;\n' +
    '}\n',
    dependencies: cxxopts_dep,
    name: 'cxxopts wrap parses --help')
  error('the cxxopts wrap did not compile a trivial --help parser')
endif

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
  'ring_census': 'test_ring_census.cpp',
  'ira_sofi': 'test_ira_sofi.cpp',
  'local_symmetry': 'test_local_symmetry.cpp',
  'classify_concurrent': 'test_classify_concurrent.cpp',
  'sphericart': 'test_sphericart.cpp',
  'cage_canon': 'test_cage_canon.cpp',
  'topo_fingerprint': 'test_topo_fingerprint.cpp',
  'ion_env': 'test_ion_env.cpp',
  'brine_front': 'test_brine_front.cpp',
  'guest_occupancy': 'test_guest_occupancy.cpp',
  'shell_rings': 'test_shell_rings.cpp',
  'ring_hop': 'test_ring_hop.cpp',
  'voronoi': 'test_voronoi.cpp',
  'structure_desc': 'test_structure_desc.cpp',
  'cage_affiliation': 'test_cage_affiliation.cpp',
  'cage_enum': 'test_cage_enum.cpp',
  'trace_threads': 'test_trace_threads.cpp',
  'pool_audit': 'test_pool_audit.cpp',
  'stencil_rows': 'test_stencil_rows.cpp',
  'ddc_regression': 'test_ddc_regression.cpp',
  'c_api': 'test_c_api.cpp',
  'brine_select': 'test_brine_select.cpp',
  'chill_offload': 'test_chill_offload.cpp',
  'chill_plus_rule': 'test_chill_plus_rule.cpp',
  'nucleation_pathway': 'test_nucleation_pathway.cpp',
  'interface_order': 'test_interface_order.cpp',
  'science_front': 'test_science_front.cpp',
  'ql_solid': 'test_ql_solid.cpp',
  'thread_scale': 'test_thread_scale.cpp',
  'thread_invariance': 'test_thread_invariance.cpp',
  'gpu_residency': 'test_gpu_residency.cpp',
  'config': 'test_config.cpp',
  'site': 'test_site.cpp',
  'wide_cell': 'test_wide_cell.cpp',
  'offset_limit': 'test_offset_limit.cpp',
  'prism_empty': 'test_prism_empty.cpp',
  'keep_axial_name': 'test_keep_axial_name.cpp',
  'basal_needle': 'test_basal_needle.cpp',
  'ddc_copies': 'test_ddc_copies.cpp',
  'hc_template': 'test_hc_template.cpp',
  'avx_gemm': 'test_avx_gemm.cpp',
  'ylm_fallback': 'test_ylm_fallback.cpp',
}

if readcon_dep.found()
  test_sources += {'readcon': 'test_readcon.cpp'}
endif
if chemfiles_dep.found()
  test_sources += {'chemfiles_io': 'test_chemfiles_io.cpp'}
  test_sources += {'chemfiles_hostile': 'test_chemfiles_hostile.cpp'}
  executable(
    'derive_confined_water_fixture',
    'fixtures/derive_confined_water_fixture.cpp',
    dependencies: chemfiles_dep,
    build_by_default: false,
  )
  if get_option('with_cli')
    test(
      'confined_water_cli',
      find_program('python3'),
      args: [
        files('test_confined_water_cli.py'),
        meson.project_build_root() / 'seams',
        meson.project_source_root() / 'input' / 'traj' / 'confined_water_40',
      ],
    )
    test(
      'molfile_cli',
      find_program('python3'),
      args: [
        files('test_molfile_cli.py'),
        meson.project_build_root() / 'seams',
        meson.project_source_root() / 'input',
      ],
    )
  endif
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
if host_machine.system() == 'darwin'
  _test_env.prepend('DYLD_LIBRARY_PATH', meson.project_build_root() / 'src')
endif
if chemfiles_libdir != ''
  _test_env.prepend('LD_LIBRARY_PATH', chemfiles_libdir)
  if host_machine.system() == 'darwin'
    _test_env.prepend('DYLD_LIBRARY_PATH', chemfiles_libdir)
  endif
endif
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
_test_rpath = (host_machine.system() == 'darwin' ? '@loader_path' : '$ORIGIN') + '/../src'
if _readcondb_rpath != ''
  _test_rpath = _test_rpath + ':' + _readcondb_rpath
endif
if chemfiles_libdir != ''
  _test_rpath = _test_rpath + ':' + chemfiles_libdir
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

if chemfiles_dep.found()
  _probe_src = meson.project_source_root() / 'cmake' / 'chemfiles_mmcif_probe.cpp'
  _configure_args = [
    files('test_chemfiles_configure.py'),
    meson.project_build_root() / 'meson-logs' / 'meson-log.txt',
    executable('chemfiles_probe_linked', _probe_src,
      cpp_args: _chemfiles_probe_args,
      dependencies: chemfiles_dep,
      link_args: _test_link_args,
      build_rpath: _test_rpath),
  ]
  if _chemfiles_system_dep.found()
    _configure_args += [executable('chemfiles_probe_system', _probe_src,
      cpp_args: _chemfiles_probe_args,
      dependencies: _chemfiles_system_dep)]
  endif
  test('chemfiles_configure', find_program('python3'),
    args: _configure_args,
    env: _test_env)
endif

# seams.h and the seams_version.h generated beside libseams.
_api_incdirs = [include_directories('../include', '../src')]

foreach name, src : test_sources
  # meson allows no nested ternary; the long-running suites get more time
  _timeout = 120
  if name == 'cage_enum' or name == 'domain' or name == 'franzblau_extra'
    # domain: a parallel run reached the last assertions after 120s.
    _timeout = 1800
  elif name == 'bulkTUM'
    # The last case is PG(2, 5). A solo coverage run is the bound.
    _timeout = 2400
  elif name == 'trace_threads'
    # GenIce cells at several widths, the mW frame, and all 148 nucleation frames.
    _timeout = 600
  elif name == 'pool_audit'
    # Ice Ic at 4096 molecules, each kernel at one thread and at eight.
    _timeout = 300
  elif name == 'cage_affiliation'
    # The seeded case passes. A parallel run was killed at 120s on the
    # last assertion.
    _timeout = 300
  elif name == 'bop' and offload_enabled
    _timeout = 600
  elif name == 'interface_order'
    _timeout = 600
  elif name == 'science_front'
    _timeout = 300
  endif
  _exe = executable('test_' + name, src,
    dependencies: [yds_dep, catch2_dep] + (name in _rang_names ? [rang_dep] : []),
    link_with: name == 'thread_invariance' ? [seams_api_lib] : [],
    cpp_args: _args,
    link_args: _test_link_args,
    include_directories: name == 'thread_invariance' ? _incdirs + _api_incdirs : _incdirs,
    build_rpath: _test_rpath)
  if name == 'thread_invariance'
    # These cases total more than the per-test budget under coverage. Keep
    # every case and width, but give each independent group the same bound.
    test(name + '_core', _exe,
      args: ['[thread_invariance]~[c_api]'],
      workdir: meson.project_source_root() / 'input',
      env: _test_env, timeout: _timeout, is_parallel: true)
    test(name + '_c_api', _exe,
      args: ['[thread_invariance][c_api]~[thread_invariance_rings]~[thread_invariance_cages]~[thread_invariance_reuse]'],
      workdir: meson.project_source_root() / 'input',
      env: _test_env, timeout: _timeout, is_parallel: true)
    test(name + '_rings', _exe,
      args: ['[thread_invariance_rings]'],
      workdir: meson.project_source_root() / 'input',
      env: _test_env, timeout: _timeout, is_parallel: true)
    test(name + '_cages', _exe,
      args: ['[thread_invariance_cages]'],
      workdir: meson.project_source_root() / 'input',
      env: _test_env, timeout: _timeout, is_parallel: true)
    test(name + '_reuse', _exe,
      args: ['[thread_invariance_reuse]'],
      workdir: meson.project_source_root() / 'input',
      env: _test_env, timeout: _timeout, is_parallel: true)
  else
    test(name, _exe,
      workdir: meson.project_source_root() / 'input',
      env: _test_env,
      timeout: _timeout,
      is_parallel: name != 'franzblau_extra' and name != 'bulkTUM',
      suite: (name == 'franzblau_extra' or name == 'bulkTUM') ? 'heavy' : [])
  endif
endforeach

test('ring_oracle',
  executable('test_ring_oracle', 'test_ring_oracle.cpp',
    dependencies: [yds_dep, catch2_dep],
    link_with: [seams_api_lib],
    cpp_args: _args,
    link_args: _test_link_args,
    include_directories: _incdirs + _api_incdirs + [include_directories('.')],
    build_rpath: _test_rpath),
  workdir: meson.project_source_root() / 'input',
  env: _test_env,
  timeout: 600)

# .lammpsdata has no chemfiles extension. The CLI, the C API read, and
# seams_path_nframes have to name the LAMMPS Data format. A malformed file
# has to be an error, not a spin in the vendored reader.
if chemfiles_dep.found() and get_option('with_cli')
  _lammpsdata_env = environment()
  _lammpsdata_env.prepend('LD_LIBRARY_PATH', meson.project_build_root() / 'src')
  if chemfiles_libdir != ''
    _lammpsdata_env.prepend('LD_LIBRARY_PATH', chemfiles_libdir)
  endif
  _lammpsdata_env.set('SEAMS_BIN', meson.project_build_root() / 'seams')
  test('lammpsdata',
    executable('test_lammpsdata', 'test_lammpsdata.cpp',
      dependencies: [catch2_dep, dependency('threads')],
      link_with: seams_api_lib,
      cpp_args: _args + ['-I' + meson.project_build_root() / 'src'],
      link_args: _test_link_args,
      include_directories: [include_directories('.'), include_directories('../include')],
      build_rpath: _test_rpath),
    env: _lammpsdata_env,
    timeout: 60)
endif

executable('interface_order', 'interface_order_main.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: _test_rpath)

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
executable('bench_ql', 'bench_ql.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: _test_rpath)

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

executable('bench_ddc_ladder', 'bench_ddc_ladder.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

executable('bench_pool_audit', 'bench_pool_audit.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

executable('bench_trace_threads', 'bench_trace_threads.cpp',
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

bench_h2h_cpc = executable('bench_h2h_cpc', 'bench_h2h_cpc.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

bench_ring_carry = executable('bench_ring_carry', 'bench_ring_carry.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

executable('bench_ring_ic', 'bench_ring_ic.cpp',
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

# Coverage arc updates are atomic. A thread pool shares those counters
# and the 1000-trial search stops scaling, so the coverage binary runs
# the same trials in child processes.
_fuzz_cpp_args = []
_fuzz_cpp_args += _args
if get_option('b_coverage')
  _fuzz_cpp_args += ['-DFUZZ_COVERAGE_PROCESSES']
endif
_fuzz_cage = executable('fuzz_cage_invariance', 'fuzz_cage_invariance.cpp',
  dependencies: [yds_dep],
  cpp_args: _fuzz_cpp_args,
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

executable('bench_hop_stride', 'bench_hop_stride.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

executable('recall_construction', 'recall_construction.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

executable('padding_probe', 'padding_probe.cpp',
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

executable('bench_rings_h2h', 'bench_rings_h2h.cpp',
  dependencies: [yds_dep],
  cpp_args: _args,
  link_args: _test_link_args,
  include_directories: _incdirs,
  build_rpath: '$ORIGIN/../src')

rings_ih_inc = [include_directories('.'), include_directories('../include')]
rings_ih_args = _args + ['-I' + meson.project_build_root() / 'src']
rings_ih_test = executable('test_rings_ih', 'test_rings_ih.cpp',
  dependencies: [catch2_dep],
  link_with: seams_api_lib,
  cpp_args: rings_ih_args,
  link_args: _test_link_args,
  include_directories: rings_ih_inc,
  build_rpath: _test_rpath)
test('rings_ih', rings_ih_test, timeout: 180, env: _test_env)

executable('bench_capi_overhead', 'bench_capi_overhead.cpp',
  dependencies: [yds_dep],
  link_with: seams_api_lib,
  cpp_args: rings_ih_args,
  link_args: _test_link_args,
  include_directories: rings_ih_inc + [_incdirs],
  build_rpath: _test_rpath)

executable('bench_rings_ih', 'bench_rings_ih.cpp',
  link_with: seams_api_lib,
  cpp_args: rings_ih_args,
  link_args: _test_link_args,
  include_directories: rings_ih_inc,
  build_rpath: _test_rpath)

test('rings_ih_external', find_program('python3'),
  args: [meson.project_source_root() / 'tests' / 'test_rings_ih_external.py'],
  timeout: 180)

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
  'mW_cutoff': ['traj/cubic-real-7.fixture.lammpstrj', '1', '1', '3.5', '2', 'cutoff'],
  'mW_knn': ['traj/cubic-real-7.fixture.lammpstrj', '1', '1', '3.5', '2', 'knn'],
  'example_cutoff': ['traj/exampleTraj.lammpstrj', '1', '2', '3.5', '2', 'cutoff'],
  'example_knn': ['traj/exampleTraj.lammpstrj', '1', '2', '3.5', '2', 'knn'],
  'levi2': ['levi:2', '2'],
  'levi3': ['levi:3', '2'],
  'levi5': ['levi:5', '1'],
}
foreach name, args : _device_host_cases
  # levi:5 is PG(2, 5). Beside the dense-frame checks it finished in
  # 5236 s. Those checks are now serial, and the limit is 6000 s.
  _host_timeout = name == 'levi5' ? 6000 : 300
  test('device_host_' + name, _check_device_host,
    args: args,
    workdir: meson.project_source_root() / 'input',
    env: _test_env,
    timeout: _host_timeout,
    is_parallel: name != 'levi5',
    suite: name == 'levi5' ? 'heavy' : [])
endforeach

# The CUDA batch kernels on the host, no device needed.
# The include is generated in gpu_host_kernels/.
test('gpu_kernels_host',
  executable('gpu_kernels_host', ['gpu_kernels_host.cpp', _cuda_kernels],
    dependencies: [yds_dep],
    cpp_args: _args,
    link_args: _test_link_args,
    include_directories: [_incdirs, _gpu_host_kernels_inc],
    build_rpath: _test_rpath),
  env: _test_env,
  timeout: 120)

_invariance_cases = {
  'mW': ['traj/cubic-real-7.fixture.lammpstrj', '1', '3.5', '1', '1'],
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
  args: ['traj/cubic-real-7.fixture.lammpstrj', '1', '3.5', '9'],
  workdir: meson.project_source_root() / 'input',
  env: _test_env,
  should_fail: true)
# 1000 translations and 1000 permutations on every ordinary frame, plus
# the eight half-box shifts. A frame above 2048 atoms is checked at one
# trial. Degree above 80 is every vertex at cap 4. Degree above 40 is
# every vertex at cap 6. The compact liquid, BCC, and cubic fixtures retain
# the dense and crystalline graph cases without requiring network access.
# The limit is 10800 s. A skip fails the run. Heavy cases
# are serial so they do not steal the Levi graph's core.
test('fuzz_cage_shift_permute', _fuzz_cage,
  args: [meson.project_source_root() / 'input' / 'traj', '1000'],
  env: _test_env,
  timeout: 10800,
  suite: 'heavy',
  is_parallel: false)
_directed_env = environment()
_directed_env.prepend('LD_LIBRARY_PATH', meson.project_build_root() / 'src')
_directed_env.set('SEAMS_LIB', seams_api_lib.full_path())
test('directed_hbonds', find_program('python3'),
  args: [files('test_directed_hbonds.py')],
  workdir: meson.project_source_root(),
  env: _directed_env,
  depends: [seams_api_lib],
  timeout: 30)
test('neighbour_record', find_program('python3'),
  args: [files('neighbour_record.py')],
  workdir: meson.project_source_root(),
  timeout: 30)
# The harness used to be a guessed build path with no dependency, so a
# flat layout or a test run that had not linked bench_h2h_cpc reported
# "missing harness or frame". The target supplies the real output path
# and is built before the test.
test('h2h_harness', find_program('python3'),
  args: [files('test_h2h_harness.py'), bench_h2h_cpc],
  workdir: meson.project_source_root(),
  env: _test_env,
  depends: [bench_h2h_cpc],
  timeout: 120)
test('h2h_campaign', find_program('python3'),
  args: [files('test_h2h_campaign.py')],
  workdir: meson.project_source_root(),
  timeout: 120)
test('snapshot_swap', find_program('python3'),
  args: [files('test_snapshot_swap.py')],
  workdir: meson.project_source_root(),
  timeout: 60)
test('bench_h2h_fetch', find_program('python3'),
  args: [files('test_bench_h2h_fetch.py')],
  workdir: meson.project_source_root(),
  timeout: 60)
test('ion_shell_rings', find_program('python3'),
  args: [files('test_ion_shell_rings.py'), meson.project_build_root() / 'seams'],
  workdir: meson.project_source_root(),
  env: _test_env,
  timeout: 120)
# Placement is sched_setaffinity, /proc/self/task and taskset: Linux only.
if host_machine.system() == 'linux'
  test('pin_host', find_program('python3'),
    args: [files('pin_host.py')],
    workdir: meson.project_source_root(),
    env: {'SEAMS_LIB': meson.project_build_root() / 'src' / 'libyodaLib.so'},
    timeout: 180,
    is_parallel: false)
endif
test('docs_sources', find_program('python3'),
  args: [files('docs_sources.py')],
  workdir: meson.project_source_root(),
  timeout: 30)
test('tracked_files_fallback', find_program('python3'),
  args: [files('test_tracked_files_fallback.py')],
  workdir: meson.project_source_root(),
  timeout: 60)
test('dist_vendor_manifest', find_program('python3'),
  args: [files('test_dist_vendor_manifest.py')],
  workdir: meson.project_source_root(),
  timeout: 60)
test('option_names', find_program('python3'),
  args: [files('option_names.py'), meson.project_build_root() / 'seams',
         meson.project_build_root()],
  workdir: meson.project_source_root(),
  timeout: 60)
test('backward_overlay', find_program('python3'),
  args: [files('backward_overlay.py')],
  workdir: meson.project_source_root(),
  timeout: 60)
test('backward_overlay_selftest', find_program('python3'),
  args: [files('backward_overlay.py'), '--self-test'],
  workdir: meson.project_source_root(),
  timeout: 60)
test('vendored_licenses', find_program('python3'),
  args: [files('vendored_licenses.py')],
  workdir: meson.project_source_root(),
  timeout: 60)
test('vendored_licenses_selftest', find_program('python3'),
  args: [files('vendored_licenses.py'), '--self-test'],
  workdir: meson.project_source_root(),
  timeout: 60)
test('cycless_pentagon', find_program('python3'),
  args: [files('test_cycless_pentagon.py'),
         seams_api_lib.full_path()],
  workdir: meson.project_source_root(),
  timeout: 120)
test('sqq_sh_triclinic', find_program('python3'),
  args: [files('test_sqq_sh.py')],
  workdir: meson.project_source_root(),
  timeout: 60)
test('ring_carry', bench_ring_carry,
  args: ['800', '2'],
  timeout: 180)
test('ovito_chill_atoms', find_program('python3'),
  args: [files('test_ovito_chill_atoms.py'),
         seams_api_lib.full_path()],
  workdir: meson.project_source_root(),
  timeout: 180)
test('ovito_ring_counts', find_program('python3'),
  args: [files('test_ovito_ring_counts.py'),
         meson.project_build_root() / 'tests' / 'bench_ring_ic'],
  workdir: meson.project_source_root(),
  timeout: 120)
test('cli_surface', find_program('python3'),
  args: [files('cli_surface.py'), meson.project_build_root() / 'seams'],
  workdir: meson.project_source_root(),
  env: {'SEAMS_HAS_CHEMFILES': chemfiles_dep.found() ? '1' : '0'},
  timeout: 600)
if chemfiles_dep.found()
  foreach _ext : ['pdb', 'gro', 'dcd']
    test('read_unreadable_' + _ext, find_program('python3'),
      args: [files('read_unreadable.py'), meson.project_build_root() / 'seams', _ext],
      workdir: meson.project_source_root(),
      timeout: 60)
  endforeach
endif
if readcon_dep.found()
  test('read_unreadable_con', find_program('python3'),
    args: [files('read_unreadable.py'), meson.project_build_root() / 'seams', 'con'],
    workdir: meson.project_source_root(),
    timeout: 60)
endif
test('mw215_deck', find_program('python3'),
  args: [files('test_mw215_deck.py')],
  workdir: meson.project_source_root(),
  timeout: 30)
test('mw215_frame', find_program('python3'),
  args: [files('test_mw215_deck.py'), 'frame'],
  workdir: meson.project_source_root(),
  timeout: 600)

test('omp_threads',
  executable('test_omp_threads', 'test_omp_threads.cpp',
    dependencies: [yds_dep, catch2_dep],
    link_with: seams_api_lib,
    cpp_args: _args + ['-I' + meson.project_build_root() / 'src'],
    link_args: _test_link_args,
    include_directories: [_incdirs, include_directories('../include')],
    build_rpath: _test_rpath),
  env: _test_env,
  timeout: 120)

# One-thread runs under taskset keep one cpu and one thread; see the driver.
if host_machine.system() == 'linux'
  test('serial_affinity', find_program('python3'),
    args: [files('test_serial_affinity.py'),
           executable('test_serial_affinity', 'test_serial_affinity.cpp',
             dependencies: [yds_dep],
             link_with: seams_api_lib,
             cpp_args: _args + ['-I' + meson.project_build_root() / 'src'],
             link_args: _test_link_args,
             include_directories: [_incdirs, include_directories('../include')],
             build_rpath: _test_rpath)],
    env: _test_env,
    timeout: 300)
endif

test('chill_mixed',
  executable('test_chill_mixed', 'test_chill_mixed.cpp',
    dependencies: [yds_dep, catch2_dep],
    link_with: seams_api_lib,
    cpp_args: _args + ['-I' + meson.project_build_root() / 'src'],
    link_args: _test_link_args,
    include_directories: [_incdirs, include_directories('../include')],
    build_rpath: _test_rpath),
  workdir: meson.project_source_root(),
  env: _test_env,
  timeout: 120)

test('trajectory',
  executable('test_trajectory', 'test_trajectory.cpp',
    dependencies: [yds_dep, catch2_dep],
    link_with: seams_api_lib,
    cpp_args: _args + ['-I' + meson.project_build_root() / 'src'],
    link_args: _test_link_args,
    include_directories: [_incdirs, include_directories('../include')],
    build_rpath: _test_rpath),
  workdir: meson.project_source_root(),
  env: _test_env,
  timeout: 120)

test('traj_capi',
  executable('test_traj_capi', 'test_traj_capi.c',
    link_with: seams_api_lib,
    include_directories: include_directories('../include'),
    c_args: ['-I' + meson.project_build_root() / 'src'],
    link_args: ['-lm'],
    build_rpath: _test_rpath),
  env: _test_env,
  timeout: 60)

test('install_pydseamslib_if_missing', find_program('python3'),
  args: [files('install_pydseamslib_if_missing.py')],
  workdir: meson.project_source_root(),
  timeout: 180)

_traj_py_env = environment()
_traj_py_env.prepend('PYTHONPATH', meson.project_source_root() / 'adapters' / 'python')
_traj_py_env.set('SEAMS_LIB', seams_api_lib.full_path())
_traj_py_env.set('SEAMS_HAS_CHEMFILES', chemfiles_dep.found() ? '1' : '0')
_traj_py_env.prepend('LD_LIBRARY_PATH', meson.project_build_root() / 'src')
if host_machine.system() == 'darwin'
  _traj_py_env.prepend('DYLD_LIBRARY_PATH', meson.project_build_root() / 'src')
endif
if chemfiles_libdir != ''
  _traj_py_env.prepend('LD_LIBRARY_PATH', chemfiles_libdir)
  if host_machine.system() == 'darwin'
    _traj_py_env.prepend('DYLD_LIBRARY_PATH', chemfiles_libdir)
  endif
endif
test('trajectory_python', find_program('python3'),
  args: ['-m', 'pytest', '-v', files('test_trajectory_dynamics.py')],
  workdir: meson.project_build_root(),
  env: _traj_py_env,
  depends: [seams_api_lib],
  timeout: 180)

# shared_module is a Mach-O bundle on Darwin, and DYLD_INSERT_LIBRARIES
# loads only a dylib. LD_PRELOAD takes the module as built.
if host_machine.system() == 'darwin'
  dump_faults = shared_library('dump_faults', 'dump_faults_preload.c',
    dependencies: dependency('dl', required: false))
else
  dump_faults = shared_module('dump_faults', 'dump_faults_preload.c', name_prefix: 'lib')
endif
dump_fault_driver = executable('dump_fault_driver', 'dump_fault_driver.cpp',
  dependencies: [yds_dep],
  include_directories: _incdirs,
  link_args: _test_link_args,
  build_rpath: _test_rpath)
_fault_wrap = find_program('dump_fault_wrap.sh')
_fault_eintr = environment()
_fault_eintr.prepend('LD_LIBRARY_PATH', meson.project_build_root() / 'src')
_fault_eintr.set('DUMP_PREAD_EINTR', '1')
_fault_eio = environment()
_fault_eio.prepend('LD_LIBRARY_PATH', meson.project_build_root() / 'src')
_fault_eio.set('DUMP_PREAD_EIO', '1')
_fault_fstat = environment()
_fault_fstat.prepend('LD_LIBRARY_PATH', meson.project_build_root() / 'src')
_fault_fstat.set('DUMP_FSTAT_EIO', '1')
test('dump_index_eintr', _fault_wrap,
  args: [dump_faults.full_path(), dump_fault_driver, 'eintr'],
  env: _fault_eintr, depends: dump_faults)
test('dump_index_eio', _fault_wrap,
  args: [dump_faults.full_path(), dump_fault_driver, 'eio'],
  env: _fault_eio, depends: dump_faults)
test('dump_fstat_keeps_atoms', _fault_wrap,
  args: [dump_faults.full_path(), dump_fault_driver, 'fstat'],
  env: _fault_fstat, depends: dump_faults)

# A consumer of libyodaLib linked without --allow-shlib-undefined and
# without sphericart on its own link line. The library the probe accepted
# has to resolve through libyodaLib plus sphericart_consumer_args, the
# -rpath-link that yds_dep also carries for a pip wheel's folder.
if sphericart_found
  _strict = cppc.has_link_argument('-Wl,--no-allow-shlib-undefined') ? ['-Wl,--no-allow-shlib-undefined'] : []
  test('sphericart_consumer',
    executable('sphericart_consumer', 'sphericart_consumer.cpp',
      link_with: ydslib,
      include_directories: include_directories('../src/include/internal'),
      link_args: _strict + sphericart_consumer_args,
      build_rpath: meson.project_build_root() / 'src'),
    timeout: 30)
endif

_wf_python = find_program('python3', required: true)
test('workflow_triggers', _wf_python,
  args: [files('check_workflow_triggers.py')],
  workdir: meson.project_source_root(),
  timeout: 30)
test('workflow_triggers_selftest', _wf_python,
  args: [files('check_workflow_triggers.py'), '--self-test'],
  workdir: meson.project_source_root(),
  timeout: 30)
