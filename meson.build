#-----------------------------------------------------------------------------------
# d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
#
# Copyright (c) 2018--present d-SEAMS core team
#
# SPDX-License-Identifier: MIT
#-----------------------------------------------------------------------------------

project('seams-core', ['cpp', 'c'],
  version: '2.0.0',
  default_options: ['warning_level=2', 'cpp_std=c++20'],
  meson_version: '>=1.3.0')

host_system = host_machine.system()

# Global accumulator variables (eOn pattern)
_args = ['-funroll-loops', '-Wno-sign-compare', '-Wno-unused-variable',
         '-Wno-unused-parameter', '-Wno-unused-but-set-variable']
_deps = []
_linkto = []
_incdirs = []

cppc = meson.get_compiler('cpp')

# GCC vs clang warning flag differences
if cppc.get_id() == 'gcc'
  _args += ['-Wno-maybe-uninitialized']
else
  _args += ['-Wno-uninitialized']
endif

# Suppress [[nodiscard]] warnings from legacy callers (to be fixed incrementally)
_args += ['-Wno-unused-result']

# Platform conditionals
if host_system == 'darwin'
  _args += ['-DOSX=TRUE']
  # The conda toolchain targets the macOS 11 SDK, whose availability
  # annotations forbid libc++ features (std::to_chars for floating point,
  # required by <format>) that the runtime libc++ we ship actually provides;
  # conda-forge's documented remedy is to disable the availability guards
  _args += ['-D_LIBCPP_DISABLE_AVAILABILITY']
endif

# Core dependencies
m_dep = cppc.find_library('m', required: false)
_deps += m_dep

eigen_dep = dependency('eigen3', version: '>=3.4.0', required: true)
_deps += [eigen_dep]

blas_dep = cppc.find_library('blas', required: true)
lapack_dep = cppc.find_library('lapack', required: true)
_deps += [blas_dep, lapack_dep]

# Shared-memory parallelism. libstdc++ routes std::execution::par through TBB,
# which conda-forge does not put in this environment, so a build relying on it
# would quietly run serially. OpenMP ships with the compiler on most toolchains.
#
# Finding the dependency is not the same as being able to link it. On the macOS
# toolchain the flags it supplies drop the default libraries from the link, and
# the failure surfaces as missing unwinder and locale symbols rather than as a
# missing libomp, so probe an actual parallel region before committing to it.
openmp_dep = dependency('openmp', required: false)
openmp_found = false
if openmp_dep.found()
  openmp_probe = '''
    #include <omp.h>
    int main() {
      int total = 0;
      #pragma omp parallel for reduction(+ : total)
      for (int i = 0; i < 64; ++i) {
        total += i;
      }
      return total == 2016 ? 0 : 1;
    }
  '''
  if cppc.links(openmp_probe, dependencies: [openmp_dep], name: 'OpenMP links')
    _deps += [openmp_dep]
    _args += ['-DSEAMS_HAS_OPENMP=1']
    openmp_found = true
  else
    message('OpenMP was found but a parallel region does not link; ' +
            'building the per-particle loops single-threaded')
  endif
endif

# MPI is optional. Analysis splits atoms across ranks and Allgathervs
# q_lm (Lechner-Dellago needs neighbours that may live on another rank).
# This is not LAMMPS-style spatial decomposition (Plimpton 1995).
# Default disabled so the CPU wheel never picks up a stray system MPI.
mpi_dep = dependency('mpi', language: 'cpp', required: get_option('with_mpi'))
mpi_found = mpi_dep.found()
if mpi_found
  _deps += [mpi_dep]
  _args += ['-DSEAMS_HAS_MPI=1']
endif

# OpenMP target offload. The probe must actually link the combined
# construct the library ships; a bare -fopenmp-targets flag or an
# empty `#pragma omp target {}` is not enough (missing nvptx bitcode
# still outlines to the host).
offload_opt = get_option('with_openmp_offload')
offload_enabled = false
if not offload_opt.disabled()
  if not openmp_found
    if offload_opt.enabled()
      error('with_openmp_offload=enabled but OpenMP does not link')
    endif
  else
    offload_arg = '-fopenmp-targets=' + get_option('offload_triple')
    offload_probe = '''
      #include <omp.h>
      int main() {
        const int n = 64;
        double x[64] = {};
        #pragma omp target teams distribute parallel for map(tofrom: x[0:n])
        for (int i = 0; i < n; ++i) x[i] = 1.0;
        return x[0] == 1.0 ? 0 : 1;
      }
    '''
    if cppc.links(offload_probe, args: ['-fopenmp', offload_arg],
                  name: 'OpenMP target offload link')
      _args += ['-DSEAMS_HAS_OFFLOAD=1', offload_arg]
      offload_enabled = true
    elif offload_opt.enabled()
      error('OpenMP target offload was requested but the compiler ' +
            'cannot link a target region for ' + get_option('offload_triple'))
    else
      message('OpenMP target offload probe failed for ' +
              get_option('offload_triple') + '; host OpenMP only')
    endif
  endif
endif

hwy_dep = dependency('hwy', version: '>=1.3.0', required: false)
if hwy_dep.found()
  _deps += [hwy_dep]
  _args += ['-DSEAMS_HAS_HWY=1']
endif

vesin_dep = dependency('vesin', required: false, fallback: ['vesin', 'vesin_dep'])
if vesin_dep.found()
  _deps += [vesin_dep]
  _args += ['-DSEAMS_HAS_VESIN=1']
endif

chemfiles_dep = dependency('chemfiles', required: false)
if chemfiles_dep.found()
  _deps += [chemfiles_dep]
  _args += ['-DSEAMS_HAS_CHEMFILES=1']
endif

readcon_dep = dependency('readcon-core', required: false,
                         fallback: ['readcon-core', 'readcon_dep'])
if readcon_dep.found()
  _deps += [readcon_dep]
  _args += ['-DSEAMS_HAS_READCON=1']
endif

# IRA/SOFI (Gunde et al.). Optional: the CPU wheel stays Horn-only
# unless the builder has libira. The C header is vendored; the probe
# must actually link libira_match.
ira_opt = get_option('with_ira')
ira_found = false
if not ira_opt.disabled()
  ira_lib_dirs = []
  # conda-forge ira installs libira.so next to ira_mod, not in $PREFIX/lib.
  python_for_ira = find_program('python3', required: false)
  if python_for_ira.found()
    ira_libdir_cmd = run_command(
      python_for_ira, '-c',
      'import pathlib, sys\n' +
      'try:\n' +
      ' import ira_mod\n' +
      ' print(pathlib.Path(ira_mod.__file__).resolve().parent / "lib")\n' +
      'except Exception:\n' +
      ' sys.exit(1)\n',
      check: false)
    if ira_libdir_cmd.returncode() == 0
      ira_lib_dirs += [ira_libdir_cmd.stdout().strip()]
    endif
  endif
  ira_dep = dependency('ira', required: false)
  if not ira_dep.found()
    ira_lib = cppc.find_library('ira', dirs: ira_lib_dirs, required: false)
    if ira_lib.found()
      ira_rpath = []
      foreach d : ira_lib_dirs
        ira_rpath += ['-Wl,-rpath,' + d]
      endforeach
      ira_dep = declare_dependency(dependencies: [ira_lib],
                                   link_args: ira_rpath)
    endif
  endif
  ira_probe = '''
    #include <iralib_interf.h>
    int main() {
      char ver[16] = {};
      long date = 0;
      libira_get_version(ver, &date);
      return 0;
    }
  '''
  ira_inc = meson.project_source_root() / 'src/include/external'
  if ira_dep.found() and cppc.links(ira_probe, dependencies: [ira_dep],
                                    args: ['-I' + ira_inc],
                                    name: 'libira C API')
    _deps += [ira_dep]
    _args += ['-DSEAMS_HAS_IRA=1']
    ira_found = true
  elif ira_opt.enabled()
    error('with_ira=enabled but libira does not link')
  endif
endif

sphericart_opt = get_option('with_sphericart')
sphericart_found = false
if not sphericart_opt.disabled()
  sph_lib_dirs = []
  sph_inc_dirs = []
  python_for_sph = find_program('python3', required: false)
  if python_for_sph.found()
    sph_dir_cmd = run_command(
      python_for_sph, '-c',
      'import pathlib, sys\n' +
      'try:\n' +
      ' import sphericart, pathlib as p\n' +
      ' root = p.Path(sphericart.__file__).resolve().parent\n' +
      ' print(root / "lib")\n' +
      ' print(root / "include")\n' +
      'except Exception:\n' +
      ' sys.exit(1)\n',
      check: false)
    if sph_dir_cmd.returncode() == 0
      lines = sph_dir_cmd.stdout().strip().split('\n')
      if lines.length() >= 2
        sph_lib_dirs += [lines[0].strip()]
        sph_inc_dirs += [lines[1].strip()]
      endif
    endif
  endif
  sphericart_dep = dependency('sphericart', required: false)
  if not sphericart_dep.found()
    sph_lib = cppc.find_library('sphericart', dirs: sph_lib_dirs, required: false)
    if sph_lib.found()
      sph_rpath = []
      foreach d : sph_lib_dirs
        sph_rpath += ['-Wl,-rpath,' + d]
      endforeach
      sph_cflags = []
      foreach d : sph_inc_dirs
        sph_cflags += ['-I' + d]
      endforeach
      sphericart_dep = declare_dependency(dependencies: [sph_lib],
                                          compile_args: sph_cflags,
                                          link_args: sph_rpath)
    endif
  endif
  sphericart_probe = '''
    #include <sphericart.hpp>
    int main() {
      sphericart::SphericalHarmonics<double> s(2);
      return 0;
    }
  '''
  sph_probe_args = []
  foreach d : sph_inc_dirs
    sph_probe_args += ['-I' + d]
  endforeach
  if sphericart_dep.found() and cppc.links(sphericart_probe,
                                          dependencies: [sphericart_dep],
                                          args: sph_probe_args,
                                          name: 'sphericart C++ API')
    _deps += [sphericart_dep]
    _args += ['-DSEAMS_HAS_SPHERICART=1']
    sphericart_found = true
  elif sphericart_opt.enabled()
    error('with_sphericart=enabled but sphericart does not link')
  endif
endif

nauty_opt = get_option('with_nauty')
nauty_found = false
if not nauty_opt.disabled()
  nauty_dep = dependency('nauty', required: false)
  if not nauty_dep.found()
    nauty_lib = cppc.find_library('nauty', required: false)
    if nauty_lib.found()
      nauty_dep = nauty_lib
    endif
  endif
  nauty_probe = '''
    #define MAXN 32
    #include <nauty/nauty.h>
    int main() {
      DEFAULTOPTIONS_GRAPH(options);
      return 0;
    }
  '''
  if nauty_dep.found() and cppc.links(nauty_probe, dependencies: [nauty_dep],
                                      name: 'nauty densenauty')
    _deps += [nauty_dep]
    _args += ['-DSEAMS_HAS_NAUTY=1']
    nauty_found = true
  elif nauty_opt.enabled()
    error('with_nauty=enabled but nauty does not link')
  endif
endif

# Report which optional backends this build actually got. Each of these
# degrades silently: without vesin the neighbour search falls back to the
# brute-force path, and without a reader the corresponding formats are simply
# absent from the library, with nothing at run time to say so.
summary({
  'SIMD distance kernel (hwy)': hwy_dep.found(),
  'Shared-memory parallelism (OpenMP)': openmp_found,
  'MPI atom split (Steinhardt)': mpi_found,
  'OpenMP target offload': offload_enabled,
  'Cell-list neighbours (vesin)': vesin_dep.found(),
  'PDB/GRO/DCD readers (chemfiles)': chemfiles_dep.found(),
  '.con reader (readcon-core)': readcon_dep.found(),
  'IRA/SOFI shape matching (libira)': ira_found,
  'Cartesian Ylm (sphericart)': sphericart_found,
  'Cage certificates (nauty)': nauty_found,
}, section: 'Optional backends', bool_yn: true)

if not readcon_dep.found()
  message('readcon-core was not found. It installs its own pkg-config file ' +
          'through cargo-c, so `cargo cinstall --prefix=$CONDA_PREFIX` in a ' +
          'readcon-core checkout is enough to enable the .con reader; a ' +
          'conda-forge recipe is in review. The subproject fallback builds ' +
          'the header with cbindgen instead and is the slower path.')
endif

# Include directories
_incdirs += include_directories([
  'src/include/internal',
  'src/include/external',
])

# Library
subdir('src')

# Tests
if get_option('with_tests')
  subdir('tests')
endif

# Python bindings
if get_option('with_python')
  subdir('python')
endif

# yodaStruct Lua CLI (optional; Lua via pkg-config or the wrap)
lua_dep = dependency('lua', required: false)
if not lua_dep.found()
  lua_dep = dependency('lua-5.4', required: false)
endif
if not lua_dep.found()
  lua_dep = dependency('lua-5.3', required: false)
endif
if not lua_dep.found() and get_option('with_lua').allowed()
  lua_dep = dependency('lua-5.4', fallback: ['lua', 'lua_dep'],
                       required: get_option('with_lua').enabled())
endif
yoda_struct = disabler()
if get_option('with_lua').allowed() and lua_dep.found()
  # Same counter-measure as the test binaries: conda LDFLAGS inject a
  # DT_RPATH to the env prefix, which outranks LD_LIBRARY_PATH and resolves a
  # stale installed libyodaLib over the freshly built one
  yoda_link_args = []
  if host_machine.system() == 'linux'
    yoda_link_args += ['-Wl,--enable-new-dtags']
  endif
  # Vendored Fennel compiler for .fnl scripts; the binary falls back to this
  # source-root path unless YODA_FENNEL_PATH overrides it at runtime.
  fennel_lua = meson.project_source_root() / 'src/include/external/fennel/fennel.lua'
  yoda_struct = executable('yodaStruct',
    files('src/main.cpp', 'src/lua_api.cpp', 'src/opt_parser.cpp'),
    dependencies: [yds_dep, lua_dep],
    cpp_args: _args + ['-DYODA_FENNEL_LUA="' + fennel_lua + '"'],
    include_directories: _incdirs,
    link_args: yoda_link_args,
    build_rpath: meson.project_build_root() / 'src',
    install: true)
  test('yoda_structure_desc', yoda_struct,
    args: ['-c', 'example_lua/structureDesc/config.yml'],
    workdir: meson.project_source_root(),
    timeout: 120)
  test('yoda_full_api', yoda_struct,
    args: ['-c', 'example_lua/full_api/config.yml'],
    workdir: meson.project_source_root(),
    timeout: 240)
  test('yoda_fennel', yoda_struct,
    args: ['-c', 'example_lua/fennel/config.yml'],
    workdir: meson.project_source_root(),
    timeout: 120)
endif

summary({
  'yodaStruct Lua CLI': lua_dep.found() and get_option('with_lua').allowed(),
}, section: 'Front ends', bool_yn: true)
