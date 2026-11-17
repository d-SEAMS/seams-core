# libseams again, with DWARF, for the abidiff step of capi_abi. -g does
# not change the generated code, so this library has the ABI of
# ../libseams.so; it keeps the same name so the sonames match. Not
# installed.
seams_abi_lib = shared_library(
  'seams',
  '../seams_api.cpp',
  dependencies: yds_dep,
  include_directories: include_directories('../../include'),
  cpp_args: _seams_api_args + ['-g'],
  link_args: _seams_api_link,
  link_depends: _seams_api_link_depends,
  build_by_default: false,
  install: false)
