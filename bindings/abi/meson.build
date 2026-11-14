# libseams_lang again, with DWARF, for capi_abi_2_13. -g does not change
# the generated code, so this library has the ABI of ../libseams_lang.so;
# it keeps the same name so the sonames match. Not installed.
seams_lang_abi_lib = shared_library(
  'seams_lang',
  files('../cxx/seams_lang.cpp'),
  dependencies: [yds_dep],
  include_directories: include_directories('../cxx'),
  cpp_args: ['-DSEAMS_LANG_VERSION="' + meson.project_version() + '"', '-g'],
  build_by_default: false,
  install: false,
)
