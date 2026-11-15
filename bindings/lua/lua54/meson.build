lua54_mod = shared_module(
  'dseams_lua54',
  files('../dseams.c'),
  name_prefix: '',
  dependencies: [lua54_dep],
  link_with: [seams_api_lib],
  include_directories: [api_inc],
  c_args: lua54_cargs,
  build_rpath: rpath,
  install: false,
)
lua54_named = custom_target(
  'dseams_core_lua54',
  input: lua54_mod,
  output: 'dseams_core.so',
  command: ['cp', '@INPUT@', '@OUTPUT@'],
  build_by_default: true,
  install: true,
  install_dir: get_option('libdir') / 'lua' / '5.4',
)
lua_mod_targets += lua54_named
lua54_prog = find_program('lua5.4', required: false)
if not lua54_prog.found()
  lua54_prog = find_program('lua', required: false)
endif
if lua54_prog.found()
  lua_runners += [lua54_prog.full_path(), meson.current_build_dir()]
endif
