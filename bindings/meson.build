# C ABI for the R and Julia fronts. Python is PydSEAMSlib and links
# libyodaLib itself. The classifiers stay in libyodaLib.

binding_oracle = executable(
  'binding_oracle',
  files('../tests/binding_oracle.cpp'),
  dependencies: [yds_dep],
  build_rpath: meson.project_build_root() / 'src',
)

seams_lang = shared_library(
  'seams_lang',
  files('cxx/seams_lang.cpp'),
  dependencies: [yds_dep],
  include_directories: include_directories('cxx'),
  cpp_args: ['-DSEAMS_LANG_VERSION="' + meson.project_version() + '"'],
  install: true,
)

seams_lang_dep = declare_dependency(
  link_with: seams_lang,
  include_directories: include_directories('cxx'),
  dependencies: [yds_dep],
)

install_headers('cxx/seams_lang.h', subdir: 'seams')

pkg = import('pkgconfig')
pkg.generate(
  name: 'seams-lang',
  filebase: 'seams-lang',
  description: 'C ABI for the d-SEAMS language fronts',
  version: meson.project_version(),
  subdirs: ['seams'],
  libraries: ['-L${libdir}', '-lseams_lang'],
  requires: ['seams-core'],
)

seams_lang_test = executable(
  'test_seams_lang',
  files('cxx/test_lang.cpp'),
  dependencies: [seams_lang_dep],
  include_directories: include_directories('cxx'),
  link_args: ['-ldl'],
)

test(
  'seams_lang_fixtures',
  seams_lang_test,
  args: [meson.project_source_root(), seams_exe.full_path()],
  depends: [seams_exe],
  timeout: 180,
)

seams_frames_test = executable(
  'test_lang_frames',
  files('cxx/test_frames.cpp'),
  dependencies: [seams_lang_dep],
  include_directories: include_directories('cxx'),
)
test('seams_lang_frames', seams_frames_test, timeout: 60)

# R is not built with the sanitizer. Loading libasan into R fails with
# "ASan runtime does not come first". The same test passes in the
# default suite, so the sanitizer run does not register it.

# Exit 77 (skipped) when pydseamslib, R dseams, or Julia is missing; CI
# passes --require-all after installing all three.
if get_option('with_cli') and get_option('with_tests') and not get_option('b_sanitize').contains('address')
  # conda-forge Julia under the toolchain loader. juliaup faults at
  # shutdown. R and pydseams are already on the default interpreter.
  chill_env = environment()
  chill_env.set('DSEAMS_ROOT', meson.project_source_root())
  _fs = import('fs')
  _julia_wrap = get_option('julia_wrap')
  if _julia_wrap != ''
    chill_env.prepend('PATH', _julia_wrap)
  else
    _julia = find_program('julia', required: false)
    if _julia.found()
      chill_env.prepend('PATH', _fs.parent(_julia.full_path()))
    endif
  endif
  _rscript = find_program('Rscript', required: false)
  if _rscript.found()
    chill_env.prepend('PATH', _fs.parent(_rscript.full_path()))
  endif
  _r_libs = get_option('r_libs_user')
  if _r_libs != ''
    chill_env.set('R_LIBS_USER', _r_libs)
  endif
  test(
    'chill_plus_fronts_agree',
    find_program('python3', required: true),
    args: [
      meson.project_source_root() / 'tests' / 'chill_plus_fronts.py',
      seams_exe.full_path(),
      seams_lang.full_path(),
    ],
    env: chill_env,
    depends: [seams_exe, seams_lang],
    timeout: 600,
  )
  test(
    'chill_plus_ctypes',
    find_program('python3', required: true),
    args: [meson.project_source_root() / 'tests' / 'test_chill_plus_ctypes.py'],
    env: {'SEAMS_LANG_LIB': seams_lang.full_path()},
    depends: [seams_lang],
    timeout: 120,
  )
  test(
    'q4_q6_chill_ctypes',
    find_program('python3', required: true),
    args: [meson.project_source_root() / 'tests' / 'test_q4_q6_chill.py'],
    env: {'SEAMS_LANG_LIB': seams_lang.full_path()},
    depends: [seams_lang],
    timeout: 120,
  )
  test(
    'wl_ice',
    find_program('python3', required: true),
    args: [meson.project_source_root() / 'tests' / 'test_wl_ice.py'],
    env: {'SEAMS_LANG_LIB': seams_lang.full_path()},
    depends: [seams_lang],
    timeout: 180,
  )
  test(
    'q4_q6_equal',
    find_program('python3', required: true),
    args: [meson.project_source_root() / 'tests' / 'test_q4_q6_equal.py'],
    env: {'SEAMS_LANG_LIB': seams_lang.full_path()},
    depends: [seams_lang],
    timeout: 180,
  )
  install_data(
    ['python/chill_plus.py', 'python/symmetry.py'],
    install_dir: get_option('datadir') / 'seams' / 'python',
  )
  # DSEAMS_PYTHON is the interpreter found at configure time. The
  # script also accepts the active conda prefix, and skips a mismatch.
  _front_py = find_program('python3', 'python', required: false)
  _front_env = environment()
  if _front_py.found()
    _front_env.set('DSEAMS_PYTHON', _front_py.full_path())
  endif
  test(
    'language_fronts',
    find_program('sh', required: true),
    args: [
      meson.project_source_root() / 'tests' / 'run_language_fronts.sh',
      meson.project_build_root(),
    ],
    env: _front_env,
    depends: [seams_lang, ydslib],
    timeout: 600,
  )
  test(
    'frames_abi',
    find_program('python3', required: true),
    args: [meson.project_source_root() / 'tests' / 'test_frames_abi.py'],
    env: {'SEAMS_LANG_LIB': seams_lang.full_path()},
    depends: [seams_lang],
    timeout: 120,
  )
  test(
    'symmetry_ctypes',
    find_program('python3', required: true),
    args: [meson.project_source_root() / 'tests' / 'test_symmetry_ctypes.py'],
    env: {'SEAMS_CORE_LIB': ydslib.full_path()},
    depends: [ydslib],
    timeout: 120,
  )
endif

# CxxWrap needs libcxxwrap_julia, which is a BinaryBuilder product and
# is not part of the pixi or Nix closure. The Julia front that ships is
# the ccall package under bindings/julia/Dseams.
subdir('lua')

cxxwrap_dep = dependency('libcxxwrap_julia', required: false)
if cxxwrap_dep.found()
  message('libcxxwrap_julia found; the shipped Julia front is still ccall')
else
  message('Julia CxxWrap is not linked. Use bindings/julia/Dseams (ccall).')
endif
