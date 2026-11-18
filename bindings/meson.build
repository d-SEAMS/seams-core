# libseams_lang: the wider C ABI the language bindings link. The
# classifiers stay in libyodaLib.

if get_option('with_tests')
  binding_oracle = executable(
    'binding_oracle',
    files('../tests/binding_oracle.cpp'),
    dependencies: [yds_dep],
    build_rpath: meson.project_build_root() / 'src',
  )

  seams_hpp_test = executable(
    'test_seams_hpp',
    files('../tests/test_seams_hpp.cpp'),
    link_with: seams_api_lib,
    include_directories: include_directories('../include'),
    cpp_args: ['-I' + meson.project_build_root() / 'src'],
  )

  test(
    'seams_hpp_match',
    seams_hpp_test,
    args: [binding_oracle.full_path(), meson.project_source_root()],
    depends: [binding_oracle, seams_api_lib],
    timeout: 180,
  )
endif

seams_lang = shared_library(
  'seams_lang',
  files('cxx/seams_lang.cpp'),
  dependencies: [yds_dep],
  include_directories: include_directories('cxx'),
  cpp_args: ['-DSEAMS_LANG_VERSION="' + meson.project_version() + '"'],
  install_rpath: '$ORIGIN',
  install: true,
)

if _abidiff.found()
  subdir('abi')
endif

seams_lang_dep = declare_dependency(
  link_with: seams_lang,
  include_directories: include_directories('cxx'),
  dependencies: [yds_dep],
)

install_headers('cxx/seams_lang.h', subdir: 'seams')

# Naming the target gives meson-uninstalled/seams-lang-uninstalled.pc the
# build-tree -L and -I. seams_lang is a shared_library(), so meson does
# not walk its dependencies into Libs.private.
pkg = import('pkgconfig')
pkg.generate(
  seams_lang,
  name: 'seams-lang',
  filebase: 'seams-lang',
  description: 'C ABI for the d-SEAMS language fronts',
  version: meson.project_version(),
  subdirs: ['seams'],
  requires: ['seams-core'],
)

if get_option('with_tests')
  seams_lang_test = executable(
    'test_seams_lang',
    files('cxx/test_lang.cpp'),
    dependencies: [seams_lang_dep],
    include_directories: include_directories('cxx'),
    link_args: ['-ldl'],
  )

  if get_option('with_cli')
    test(
      'seams_lang_fixtures',
      seams_lang_test,
      args: [meson.project_source_root(), seams_exe.full_path()],
      depends: [seams_exe],
      timeout: 180,
    )
  endif

  seams_frames_test = executable(
    'test_lang_frames',
    files('cxx/test_frames.cpp'),
    dependencies: [seams_lang_dep],
    include_directories: include_directories('cxx'),
  )
  test('seams_lang_frames', seams_frames_test, timeout: 60)

  # ctypes checks of libseams_lang and libyodaLib. The language bindings
  # live in their own repositories; see MAINTAINERS/repositories.org.
  if get_option('with_cli') and not get_option('b_sanitize').contains('address')
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
    _tip4p_env = environment()
    _tip4p_env.set('SEAMS_LANG_LIB', seams_lang.full_path())
    _tip4p_env.prepend('LD_LIBRARY_PATH', meson.project_build_root() / 'src')
    if chemfiles_libdir != ''
      _tip4p_env.prepend('LD_LIBRARY_PATH', chemfiles_libdir)
    endif
    test(
      'tip4p_overhead',
      find_program('python3', required: true),
      args: [meson.project_source_root() / 'tests' / 'test_tip4p_overhead.py'],
      workdir: meson.project_source_root(),
      env: _tip4p_env,
      depends: [seams_lang],
      timeout: 180,
    )
  endif
endif
