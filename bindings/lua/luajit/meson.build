luajit_mod = shared_module(
  'dseams_luajit',
  files('../dseams.c'),
  name_prefix: '',
  dependencies: [luajit_dep],
  link_with: [seams_api_lib],
  include_directories: [api_inc],
  c_args: [version_inc],
  build_rpath: rpath,
  install: false,
)
luajit_named = custom_target(
  'dseams_core_luajit',
  input: luajit_mod,
  output: 'dseams_core.so',
  command: ['cp', '@INPUT@', '@OUTPUT@'],
  build_by_default: true,
  install: true,
  install_dir: get_option('libdir') / 'lua' / '5.1',
)
lua_mod_targets += luajit_named
lua_runners += [find_program('luajit').full_path(), meson.current_build_dir()]
