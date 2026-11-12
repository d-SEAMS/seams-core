# Engine drivers. They call the C API and link the distribution
# LAMMPS and GROMACS, which the module compiler cannot link.

link_host = find_program('link_host.sh')
fs = import('fs')
seams_libdir = fs.parent(seams_api_lib.full_path())

run_hook = custom_target(
  'run_hook',
  input: files('lammps/run_hook.cpp'),
  output: 'run_hook',
  command: [
    link_host, '-O2', '-std=c++17',
    '-I/usr/include/lammps',
    '-I' + meson.project_source_root() / 'include',
    '-I' + meson.project_build_root() / 'src',
    '-o', '@OUTPUT@', '@INPUT@',
    '-L/usr/lib/x86_64-linux-gnu',
    '-Wl,--allow-shlib-undefined',
    '-Wl,-rpath,/usr/lib/x86_64-linux-gnu',
    '-Wl,-rpath,' + seams_libdir,
    '-llammps', '/usr/lib/x86_64-linux-gnu/libcurl.so.4',
    seams_api_lib.full_path(),
  ],
  depends: [seams_api_lib],
  build_by_default: true,
)

test(
  'lammps_hook',
  find_program('env'),
  args: [
    '-u', 'LIBRARY_PATH',
    'LD_LIBRARY_PATH=' + seams_libdir,
    find_program('sh').full_path(),
    '-c', 'cd "$1" && exec "$2" in.dseams', 'lammps',
    meson.project_source_root() / 'engines/lammps',
    run_hook.full_path(),
  ],
  depends: [run_hook],
  timeout: 120,
)

seams_ice = custom_target(
  'seams_ice',
  input: files('gromacs/seams_ice.cpp'),
  output: 'seams_ice',
  command: [
    link_host, '-O2', '-std=c++17',
    '-I/usr/include', '-DGMX_DOUBLE=0', '-fopenmp',
    '-I' + meson.project_source_root() / 'include',
    '-I' + meson.project_build_root() / 'src',
    '-o', '@OUTPUT@', '@INPUT@',
    '-L/usr/lib/x86_64-linux-gnu',
    '-Wl,--allow-shlib-undefined',
    '-Wl,-rpath,/usr/lib/x86_64-linux-gnu',
    '-Wl,-rpath,' + seams_libdir,
    '-lgromacs', '-lm', '-lfftw3f',
    seams_api_lib.full_path(),
  ],
  depends: [seams_api_lib],
  build_by_default: true,
)

test(
  'gromacs_mdrun',
  find_program('gromacs/run_mdrun.sh'),
  args: [seams_ice.full_path(), seams_libdir],
  depends: [seams_ice],
  timeout: 180,
)
