# Standalone timing driver. Not part of the Catch2 suite.
bench_link = _link_args
if host_machine.system() == 'linux'
  bench_link += ['-Wl,--enable-new-dtags']
endif

seams_bench = executable(
  'seams_bench',
  files('analyse.cpp'),
  dependencies: [yds_dep],
  include_directories: _incdirs,
  cpp_args: _args,
  link_args: bench_link,
  build_rpath: meson.project_build_root() / 'src',
  install: false,
)
