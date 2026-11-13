# Fortran 2008 module. fpm.toml in this directory is the other build.

fortran_opt = get_option('with_fortran')
if fortran_opt.disabled()
  subdir_done()
endif

fc_added = add_languages('fortran', required: fortran_opt, native: false)
if not fc_added
  subdir_done()
endif

fortran_lib_src = files(
  'src/seams_kinds.f90',
  'src/seams_c.f90',
  'src/seams_support.f90',
  'src/seams.f90',
  'src/seams_methods.f90',
)

if get_option('with_tests')
  fortran_env = environment()
  fortran_env.set('OMP_NUM_THREADS', '1')
  capi = executable(
    'test_fortran_capi',
    files('test_seams.f90', 'src/seams_c.f90'),
    link_with: seams_api_lib,
  )
  test(
    'fortran_capi',
    capi,
    args: [meson.project_source_root()],
    env: fortran_env,
    timeout: 180,
  )
  high = executable(
    'test_fortran_high',
    files('test/test_high.f90') + fortran_lib_src,
    link_with: seams_api_lib,
  )
  test(
    'fortran_high',
    high,
    args: [meson.project_source_root()],
    env: fortran_env,
    timeout: 180,
  )
  cubic = executable(
    'fortran_cubic',
    files('example/cubic.f90') + fortran_lib_src,
    link_with: seams_api_lib,
  )
  test('fortran_cubic', cubic, env: fortran_env, timeout: 180)
endif
