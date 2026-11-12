# C ABI for the R and Julia fronts. Python is PydSEAMSlib and links
# libyodaLib itself. The classifiers stay in libyodaLib.

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

# Exit 77 (skipped) when pydseamslib, R dseams, or Julia is missing; CI
# passes --require-all after installing all three.
if get_option('with_cli') and get_option('with_tests')
  test(
    'chill_plus_fronts_agree',
    find_program('python3', required: true),
    args: [
      meson.project_source_root() / 'tests' / 'chill_plus_fronts.py',
      seams_exe.full_path(),
      seams_lang.full_path(),
    ],
    env: {'DSEAMS_ROOT': meson.project_source_root()},
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
  install_data(
    'python/chill_plus.py',
    install_dir: get_option('datadir') / 'seams' / 'python',
  )
endif

# CxxWrap needs libcxxwrap_julia, which is a BinaryBuilder product and
# is not part of the pixi or Nix closure. The Julia front that ships is
# the ccall package under bindings/julia/Dseams.
cxxwrap_dep = dependency('libcxxwrap_julia', required: false)
if cxxwrap_dep.found()
  message('libcxxwrap_julia found; the shipped Julia front is still ccall')
else
  message('Julia CxxWrap is not linked. Use bindings/julia/Dseams (ccall).')
endif
