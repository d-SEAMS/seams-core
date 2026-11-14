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

# seams_rings_find per frame through the public C API (bench/ring_frames.cpp).
seams_ring_frames = executable(
  'seams_ring_frames',
  files('ring_frames.cpp'),
  link_with: seams_api_lib,
  include_directories: include_directories('../include'),
  cpp_args: ['-I' + meson.project_build_root() / 'src'],
  install: false,
)
