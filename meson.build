#-----------------------------------------------------------------------------------
# d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
#
# Copyright (c) 2018--present d-SEAMS core team
#
# SPDX-License-Identifier: MIT
#-----------------------------------------------------------------------------------

project('seams-core', ['cpp', 'c'],
  version: '2.9.0',
  default_options: ['warning_level=2', 'cpp_std=c++20'],
  meson_version: '>=1.3.0')

host_system = host_machine.system()

# Global accumulator variables (eOn pattern)
_args = ['-funroll-loops']
_link_args = []
_deps = []
_linkto = []
_incdirs = []

cppc = meson.get_compiler('cpp')
cpp_id = cppc.get_id()

# GCC/clang warning suppressions. nvc++ (id nvidia/pgi) rejects several
# of these -Wno flags, so they stay off that toolchain.
if cpp_id == 'gcc' or cpp_id == 'clang'
  _args += ['-Wno-sign-compare', '-Wno-unused-variable',
            '-Wno-unused-parameter', '-Wno-unused-but-set-variable']
  if cpp_id == 'gcc'
    _args += ['-Wno-maybe-uninitialized']
  else
    _args += ['-Wno-uninitialized']
  endif
  # Suppress [[nodiscard]] warnings from legacy callers
  _args += ['-Wno-unused-result']
endif

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

# Prefer FlexiBLAS (one wrapper for BLAS and LAPACK). The conda-forge
# ABI packages still provide libblas/liblapack; those are the fallback
# when pkg-config does not find FlexiBLAS.
flexiblas_dep = dependency('flexiblas', required: false)
if not flexiblas_dep.found()
  flexiblas_dep = cppc.find_library('flexiblas', required: false)
endif
if flexiblas_dep.found()
  _deps += flexiblas_dep
  blas_backend = 'flexiblas'
else
  blas_dep = cppc.find_library('blas', required: true)
  lapack_dep = cppc.find_library('lapack', required: true)
  _deps += [blas_dep, lapack_dep]
  blas_backend = 'blas+lapack'
endif

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
# still outlines to the host). nvc++ takes -mp=gpu; clang/LLVM takes
# -fopenmp-targets=<triple>. The first set that links wins, and those
# flags go on both the compile and the link line.
offload_opt = get_option('with_openmp_offload')
offload_enabled = false
offload_flags = []
if not offload_opt.disabled()
  if not openmp_found
    if offload_opt.enabled()
      error('with_openmp_offload=enabled but OpenMP does not link')
    endif
  else
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
    offload_triple = get_option('offload_triple')
    offload_try = [
      ['-mp=gpu'],
      ['-mp=gpu', '-gpu=cc80'],
      ['-fopenmp', '-fopenmp-targets=' + offload_triple],
    ]
    foreach flags : offload_try
      if not offload_enabled
        probe_name = 'OpenMP target offload'
        foreach f : flags
          probe_name += ' ' + f
        endforeach
        if cppc.links(offload_probe, args: flags, dependencies: [openmp_dep],
                      name: probe_name)
          offload_flags = flags
          offload_enabled = true
        endif
      endif
    endforeach
    if offload_enabled
      _args += ['-DSEAMS_HAS_OFFLOAD=1']
      _args += offload_flags
      _link_args += offload_flags
    elif offload_opt.enabled()
      error('OpenMP target offload was requested but the compiler ' +
            'cannot link a target region (tried -mp=gpu and -fopenmp-targets=' +
            offload_triple + ')')
    else
      message('OpenMP target offload probe failed for -mp=gpu and ' +
              offload_triple + '; host OpenMP only')
    endif
  endif
endif

hwy_dep = dependency('hwy', version: '>=1.3.0', required: false)
if hwy_dep.found()
  _deps += [hwy_dep]
  _args += ['-DSEAMS_HAS_HWY=1']
endif

gpulite_opt = get_option('with_gpulite')
gpulite_dep = dependency('gpulite', required: false,
                         fallback: ['gpulite', 'gpulite_dep'])
if gpulite_opt.enabled() and not gpulite_dep.found()
  error('with_gpulite=enabled but gpulite was not found')
endif
if gpulite_dep.found() and not gpulite_opt.disabled()
  _deps += [gpulite_dep]
  _args += ['-DSEAMS_HAS_GPULITE=1']
endif

vesin_dep = dependency('vesin', required: false, fallback: ['vesin', 'vesin_dep'])
if vesin_dep.found()
  _deps += [vesin_dep]
  _args += ['-DSEAMS_HAS_VESIN=1']
endif

_linkcell_opts = []
if gpulite_dep.found() and not gpulite_opt.disabled()
  _linkcell_opts += ['with_gpulite=enabled']
endif
linkcell_dep = dependency('linkcell', required: false,
                          fallback: ['linkcell', 'linkcell_dep'],
                          default_options: _linkcell_opts)
if linkcell_dep.found()
  _deps += [linkcell_dep]
  _args += ['-DSEAMS_HAS_LINKCELL=1']
endif

chemfiles_dep = dependency('chemfiles', required: false)
if chemfiles_dep.found()
  _deps += [chemfiles_dep]
  _args += ['-DSEAMS_HAS_CHEMFILES=1']
endif

# The Rust reader follows the engine's library kind, so a static engine
# embeds it and a shared engine ships it alongside; a shared reader under a
# static engine would point consumers at a build-directory path.
_readcon_kind = get_option('default_library') == 'static' ? 'static' : 'shared'
readcon_dep = dependency('readcon-core', required: false,
                         fallback: ['readcon-core', 'readcon_dep'],
                         default_options: ['default_library=' + _readcon_kind])
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
# degrades silently: without vesin the cutoff neighbour search falls
# back to the brute-force path, without linkcell the k-nearest walk
# uses the in-tree cell list, and without a reader the corresponding
# formats are simply absent from the library, with nothing at run time
# to say so.
summary({
  'BLAS/LAPACK': blas_backend,
  'SIMD distance kernel (hwy)': hwy_dep.found(),
  'Shared-memory parallelism (OpenMP)': openmp_found,
  'MPI atom split (Steinhardt)': mpi_found,
  'OpenMP target offload': offload_enabled,
  'Device-resident batches (gpulite)': gpulite_dep.found() and not gpulite_opt.disabled(),
  'Cell-list neighbours (vesin)': vesin_dep.found(),
  'Periodic k-nearest (linkcell)': linkcell_dep.found(),
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

# Python bindings are not in this repository.
if get_option('with_python')
  error('Python is pydseams: https://github.com/d-SEAMS/PydSEAMSlib')
endif

# Lua is a library in the yodaStruct / luadseams repo, not this one.
if not get_option('with_lua').disabled()
  error('Lua/Fennel is the dseams library: https://github.com/d-SEAMS/yodaStruct')
endif

if get_option('with_cli')
  seams_cli_link = _link_args
  if host_machine.system() == 'linux'
    seams_cli_link += ['-Wl,--enable-new-dtags']
  endif
  seams_exe = executable(
    'seams',
    files('src/seams_cli.cpp'),
    dependencies: [yds_dep],
    include_directories: _incdirs,
    cpp_args: _args + ['-DSEAMS_VERSION="' + meson.project_version() + '"'],
    link_args: seams_cli_link,
    install: true)
  test(
    'seams_read',
    seams_exe,
    args: ['read', 'input/traj/exampleTraj.lammpstrj'],
    workdir: meson.project_source_root(),
    timeout: 60)
  test('seams_help', seams_exe, args: ['--help'])
  test('seams_version', seams_exe, args: ['--version'])
  test('seams_features', seams_exe, args: ['--features'])
  test('seams_print_config', seams_exe, args: ['--print-config'])

  python3 = find_program('python3', required: true)
  test(
    'seams_frame_range_preserves_order',
    python3,
    args: [
      '-c',
      'import re,subprocess,sys; p=subprocess.run([sys.argv[1],"--frame","1","--last","3","--jobs","4","read",sys.argv[2]],capture_output=True,text=True); frames=[int(m.group(1)) for m in re.finditer(r"frame (\\d+)",p.stdout)]; sys.exit(0 if p.returncode==0 and frames==[1,2,3] else 1)',
      seams_exe.full_path(),
      meson.project_source_root() / 'input' / 'traj' / 'mW_cubic.lammpstrj',
    ],
    depends: [seams_exe],
    timeout: 60)
  test(
    'seams_single_frame_format_range_is_bounded',
    python3,
    args: [
      '-c',
      'import re,subprocess,sys; p=subprocess.run([sys.argv[1],"--frame","1","--last","2","read",sys.argv[2]],capture_output=True,text=True); frames=[int(m.group(1)) for m in re.finditer(r"frame (\\d+)",p.stdout)]; sys.exit(0 if p.returncode==0 and frames==[1] else 1)',
      seams_exe.full_path(),
      meson.project_source_root() / 'input' / 'xyz' / 'configuration.xyz',
    ],
    depends: [seams_exe],
    timeout: 60)
  test(
    'seams_json_output_is_versioned',
    python3,
    args: [
      '-c',
      'import json,subprocess,sys; p=subprocess.run([sys.argv[1],"--format","json","read",sys.argv[2]],capture_output=True,text=True); r=json.loads(p.stdout); sys.exit(0 if p.returncode==0 and r["schema"]=="dseams.cli/v1" and r["command"]=="read" and r["frame"]==1 and r["status"]==0 and r["text"].startswith("nop ") else 1)',
      seams_exe.full_path(),
      meson.project_source_root() / 'input' / 'traj' / 'exampleTraj.lammpstrj',
    ],
    depends: [seams_exe],
    timeout: 60)
  test(
    'seams_hbonds_reports_selected_atoms',
    python3,
    args: [
      '-c',
      'import json,subprocess,sys; p=subprocess.run([sys.argv[1],"--format","json","hbonds",sys.argv[2],"--type","2","--htype","1"],capture_output=True,text=True); r=json.loads(p.stdout); sys.exit(0 if p.returncode==0 and r["text"].startswith("nop 250 hbonds ") else 1)',
      seams_exe.full_path(),
      meson.project_source_root() / 'input' / 'traj' / 'exampleTraj.lammpstrj',
    ],
    depends: [seams_exe],
    timeout: 60)
  test(
    'seams_json_range_is_ordered',
    python3,
    args: [
      '-c',
      'import json,subprocess,sys; p=subprocess.run([sys.argv[1],"--format","json","--frame","1","--last","3","--jobs","4","read",sys.argv[2]],capture_output=True,text=True); rows=[json.loads(line) for line in p.stdout.splitlines() if line]; sys.exit(0 if p.returncode==0 and [r["frame"] for r in rows]==[1,2,3] and all(r["schema"]=="dseams.cli/v1" for r in rows) else 1)',
      seams_exe.full_path(),
      meson.project_source_root() / 'input' / 'traj' / 'mW_cubic.lammpstrj',
    ],
    depends: [seams_exe],
    timeout: 60)
  test(
    'seams_strict_input_rejects_missing_file',
    python3,
    args: [
      '-c',
      'import subprocess,sys; p=subprocess.run([sys.argv[1],"--strict-input","read","/definitely/missing/seams-input.lammpstrj"],capture_output=True,text=True); sys.exit(0 if p.returncode==2 and "contains no atoms" in p.stdout else 1)',
      seams_exe.full_path(),
    ],
    depends: [seams_exe],
    timeout: 60)
  test(
    'seams_help_mentions_rdf',
    python3,
    args: [
      '-c',
      'import subprocess,sys; p=subprocess.run([sys.argv[1],"--help"],capture_output=True,text=True); sys.exit(0 if p.returncode==0 and "rdf" in p.stdout+p.stderr else 1)',
      seams_exe.full_path(),
    ],
    depends: [seams_exe])
  test(
    'seams_help_mentions_density_z',
    python3,
    args: [
      '-c',
      'import subprocess,sys; p=subprocess.run([sys.argv[1],"--help"],capture_output=True,text=True); text=p.stdout+p.stderr; sys.exit(0 if p.returncode==0 and "density-z" in text and "--axis" in text else 1)',
      seams_exe.full_path(),
    ],
    depends: [seams_exe])
  test(
    'seams_density_z',
    python3,
    args: [
      '-c',
      'import os,subprocess,sys; env=os.environ.copy(); lib=os.path.join(os.path.dirname(sys.argv[1]),"src"); env["LD_LIBRARY_PATH"]=lib+(":"+env["LD_LIBRARY_PATH"] if env.get("LD_LIBRARY_PATH") else ""); p=subprocess.run([sys.argv[1],"density-z",sys.argv[2],"--bins","4","--type","0"],capture_output=True,text=True,env=env); lines=[ln for ln in p.stdout.splitlines() if ln.strip()]; sys.exit(0 if p.returncode==0 and lines and lines[0]=="# z rho" and len(lines)==5 else 1)',
      seams_exe.full_path(),
      meson.project_source_root() / 'input' / 'traj' / 'exampleTraj.lammpstrj',
    ],
    depends: [seams_exe])
  test(
    'seams_unknown_command_exit_2',
    python3,
    args: [
      '-c',
      'import subprocess,sys; p=subprocess.run([sys.argv[1],"nope","x"],capture_output=True,text=True); sys.exit(0 if p.returncode==2 else 1)',
      seams_exe.full_path(),
    ],
    depends: [seams_exe])
  test(
    'seams_chill_family_ionicLiquid_exit_2',
    python3,
    args: [
      '-c',
      'import subprocess,sys; p=subprocess.run([sys.argv[1],"chill",sys.argv[2],"--family","ionicLiquid"],capture_output=True,text=True); text=p.stderr+p.stdout; sys.exit(0 if p.returncode==2 and "ionicLiquid" in text else 1)',
      seams_exe.full_path(),
      meson.project_source_root() / 'input' / 'traj' / 'exampleTraj.lammpstrj',
    ],
    depends: [seams_exe])
  test(
    'seams_cages_family_des_exit_2',
    python3,
    args: [
      '-c',
      'import subprocess,sys; p=subprocess.run([sys.argv[1],"cages",sys.argv[2],"--family","des"],capture_output=True,text=True); text=p.stderr+p.stdout; sys.exit(0 if p.returncode==2 and "des" in text else 1)',
      seams_exe.full_path(),
      meson.project_source_root() / 'input' / 'traj' / 'exampleTraj.lammpstrj',
    ],
    depends: [seams_exe])
  test(
    'seams_chill_exampleTraj',
    seams_exe,
    args: ['chill', 'input/traj/exampleTraj.lammpstrj'],
    workdir: meson.project_source_root(),
    timeout: 120)
  test(
    'seams_help_mentions_domains',
    python3,
    args: [
      '-c',
      'import subprocess,sys; p=subprocess.run([sys.argv[1],"--help"],capture_output=True,text=True); text=p.stdout+p.stderr; sys.exit(0 if p.returncode==0 and "domains" in text and "--subset" in text else 1)',
      seams_exe.full_path(),
    ],
    depends: [seams_exe])
  test(
    'seams_domains_needs_site',
    python3,
    args: [
      '-c',
      'import subprocess,sys; p=subprocess.run([sys.argv[1],"domains",sys.argv[2]],capture_output=True,text=True); sys.exit(0 if p.returncode==2 else 1)',
      seams_exe.full_path(),
      meson.project_source_root() / 'input' / 'traj' / 'exampleTraj.lammpstrj',
    ],
    depends: [seams_exe])
  test(
    'seams_help_mentions_family',
    python3,
    args: [
      '-c',
      'import subprocess,sys; p=subprocess.run([sys.argv[1],"--help"],capture_output=True,text=True); text=p.stdout+p.stderr; sys.exit(0 if p.returncode==0 and "--family" in text and "waterIce" in text else 1)',
      seams_exe.full_path(),
    ],
    depends: [seams_exe])
  test(
    'seams_cn_site_site',
    python3,
    args: [
      '-c',
      'import subprocess,sys; p=subprocess.run([sys.argv[1],"cn",sys.argv[2],"--types","2,2","--cutoff","3.5"],capture_output=True,text=True); text=p.stdout+p.stderr; sys.exit(0 if p.returncode==0 and "site-site" in text else 1)',
      seams_exe.full_path(),
      meson.project_source_root() / 'input' / 'traj' / 'exampleTraj.lammpstrj',
    ],
    depends: [seams_exe],
    timeout: 60)
  test(
    'seams_cn_ions_cage',
    python3,
    args: [
      '-c',
      'import subprocess,sys; p=subprocess.run([sys.argv[1],"cn",sys.argv[2],"--ions","--site","1=cationHead,2=anion","--cutoff","3.5"],capture_output=True,text=True); text=p.stdout+p.stderr; sys.exit(0 if p.returncode==0 and "cage" in text else 1)',
      seams_exe.full_path(),
      meson.project_source_root() / 'input' / 'traj' / 'exampleTraj.lammpstrj',
    ],
    depends: [seams_exe],
    timeout: 60)
  test(
    'seams_cages_signature_guests',
    python3,
    args: [
      '-c',
      'import subprocess,sys; p=subprocess.run([sys.argv[1],"cages",sys.argv[2],"--graph","cutoff","--cutoff","3.5","--type","1","--signature","512","--guest-types","1","--guest-radius","4.0"],capture_output=True,text=True); text=p.stdout+p.stderr; sys.exit(0 if p.returncode==0 and "cages" in text and "guests" in text and "occupied" in text else 1)',
      seams_exe.full_path(),
      meson.project_source_root() / 'input' / 'traj' / 'genice_sI.lammpstrj',
    ],
    depends: [seams_exe],
    timeout: 120)
  test(
    'seams_cages_per_atom_dump',
    python3,
    args: [
      '-c',
      'import subprocess,sys,tempfile,os; d=tempfile.mkdtemp(); out=os.path.join(d,"cages.lammpstrj"); p=subprocess.run([sys.argv[1],"cages",sys.argv[2],"--graph","cutoff","--cutoff","3.5","--type","1","--per-atom",out],capture_output=True,text=True); text=open(out).read() if os.path.exists(out) else ""; lines=text.splitlines(); n=int(lines[3]) if len(lines)>3 else -1; rows=[l for l in lines[9:9+n]]; sys.exit(0 if p.returncode==0 and "ITEM: ATOMS id type x y z cage" in text and len(rows)==n and all(l.split()[-1] in ("0","1","2","3") for l in rows) else 1)',
      seams_exe.full_path(),
      meson.project_source_root() / 'input' / 'traj' / 'genice_sI.lammpstrj',
    ],
    depends: [seams_exe],
    timeout: 120)
  test(
    'seams_pairs_contact',
    python3,
    args: [
      '-c',
      'import subprocess,sys; p=subprocess.run([sys.argv[1],"pairs",sys.argv[2],"--site","1=cationHead,2=anion"],capture_output=True,text=True); text=p.stdout+p.stderr; sys.exit(0 if p.returncode==0 and "contact-pair" in text else 1)',
      seams_exe.full_path(),
      meson.project_source_root() / 'input' / 'traj' / 'exampleTraj.lammpstrj',
    ],
    depends: [seams_exe],
    timeout: 60)
endif
