#-----------------------------------------------------------------------------------
# d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
#
# Copyright (c) 2018--present d-SEAMS core team
#
# SPDX-License-Identifier: MIT
#-----------------------------------------------------------------------------------

project('seams-core', ['cpp', 'c'],
  version: '2.14.0',
  default_options: ['warning_level=2', 'cpp_std=c++20', 'b_ndebug=if-release'],
  meson_version: '>=1.3.0')

_capi_san = get_option('capi_sanitize')
if _capi_san != 'none'
  _san_args = ['-fsanitize=' + _capi_san, '-fno-omit-frame-pointer']
  if _capi_san.contains('undefined')
    # GCC leaves float-to-int overflow out of -fsanitize=undefined.
    _san_args += ['-fsanitize=float-cast-overflow', '-fno-sanitize-recover=all']
  endif
  add_project_arguments(_san_args, language: ['c', 'cpp'])
  add_project_link_arguments(_san_args, language: ['c', 'cpp'])
endif

# libFuzzer coverage for this project's code only. Subprojects stay
# uninstrumented, since their examples link no sanitizer runtime.
if get_option('capi_fuzz_coverage')
  if meson.get_compiler('c').get_id() != 'clang' or meson.get_compiler('cpp').get_id() != 'clang'
    error('capi_fuzz_coverage needs clang for C and C++')
  endif
  add_project_arguments('-fsanitize=fuzzer-no-link', language: ['c', 'cpp'])
endif
# Clang links a sanitizer runtime into executables only, so an
# instrumented shared library leaves the runtime's symbols for the
# executable to bind. -z undefs comes after b_lundef's --no-undefined and
# allows that.
if meson.get_compiler('cpp').get_id() == 'clang' and (_capi_san != 'none' or get_option('capi_fuzz_coverage'))
  add_project_link_arguments('-Wl,-z,undefs', language: ['c', 'cpp'])
endif

host_system = host_machine.system()

_dist_python = find_program('python3', 'python', required: false)
if _dist_python.found()
  meson.add_dist_script(_dist_python, files('scripts/dist_strip_vcs.py'))
endif

# The wrap overlay's backward.hpp must stay the header this tree
# vendored. Meson skips wrap patches when the subproject directory
# already has a meson.build, which is how the offline recipes pre-seed
# it, so this check runs at configure rather than only inside a test.
_carry_python = find_program('python3', native: true, required: true)
run_command(
  _carry_python,
  files('tests/backward_overlay.py'),
  check: true,
)

# Global accumulator variables (eOn pattern)
_args = []
_link_args = []
_deps = []
_linkto = []
_incdirs = []

cppc = meson.get_compiler('cpp')
cpp_id = cppc.get_id()

# GCC/clang warning suppressions and -funroll-loops. nvc++ (id
# nvidia / nvidia_hpc / pgi) rejects those switches.
if cpp_id == 'gcc' or cpp_id == 'clang'
  _args += ['-funroll-loops']
  _args += ['-Wno-sign-compare', '-Wno-unused-variable',
            '-Wno-unused-parameter', '-Wno-unused-but-set-variable']
  if cpp_id == 'gcc'
    _args += ['-Wno-maybe-uninitialized']
  else
    _args += ['-Wno-uninitialized']
    # Signed index conversions, shadowed names, and unknown attributes
    # are errors on Clang. The flags stay on this project; subprojects
    # keep their own warning set.
    _args += ['-Wsign-conversion', '-Wshorten-64-to-32', '-Wshadow',
              '-Wunknown-attributes', '-Werror=sign-conversion',
              '-Werror=shorten-64-to-32', '-Werror=shadow',
              '-Werror=unknown-attributes']
  endif
  # Suppress [[nodiscard]] warnings from legacy callers
  _args += ['-Wno-unused-result']
endif

# A non-atomic arc update races when a test uses more than one thread.
# The counter goes negative, and the report then drops the line.
if get_option('b_coverage') and (cpp_id == 'gcc' or cpp_id == 'clang')
  _args += ['-fprofile-update=atomic']
  # Subprojects also have threaded kernels, so their counters must not
  # use a read-modify-write update shared by several workers.
  add_global_arguments('-fprofile-update=atomic', language: ['c', 'cpp'])
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

# The wheel build falls back to the Eigen subproject. Its headers are
# not this project's sources, so -Werror=sign-conversion stays on our
# code and does not fire inside Eigen.
eigen_dep = dependency('eigen3', version: '>=3.4.0', include_type: 'system',
                       required: true)
_deps += [eigen_dep]

# Prefer FlexiBLAS (one wrapper for BLAS and LAPACK). The conda-forge
# ABI packages still provide libblas/liblapack; those are the fallback
# when pkg-config does not find FlexiBLAS.
# No source calls BLAS or LAPACK today; disabling them keeps the wheel
# from grafting a BLAS it never uses.
blas_opt = get_option('with_blas')
blas_backend = 'none'
flexiblas_dep = dependency('flexiblas', required: blas_opt.disabled() ? blas_opt : false)
if not flexiblas_dep.found() and not blas_opt.disabled()
  flexiblas_dep = cppc.find_library('flexiblas', required: false)
endif
if flexiblas_dep.found()
  _deps += flexiblas_dep
  blas_backend = 'flexiblas'
  # libflexiblas.so.3 is not on the default loader path. The link
  # line is shared by every binary, so the directory is recorded here.
  _fb_file = run_command(cppc, '-print-file-name=libflexiblas.so', check: false)
  if _fb_file.returncode() == 0
    _fb_path = _fb_file.stdout().strip()
    if _fb_path != '' and _fb_path != 'libflexiblas.so'
      _fb_dir = import('fs').parent(_fb_path)
      _link_args += ['-Wl,-rpath,' + _fb_dir]
      # -rpath-link is a GNU ld option. Apple ld64 rejects it.
      if host_machine.system() != 'darwin'
        _link_args += ['-Wl,-rpath-link,' + _fb_dir]
      endif
    endif
  endif
elif not blas_opt.disabled()
  blas_dep = cppc.find_library('blas', required: blas_opt)
  lapack_dep = cppc.find_library('lapack', required: blas_opt)
  if blas_dep.found() and lapack_dep.found()
    _deps += [blas_dep, lapack_dep]
    blas_backend = 'blas+lapack'
  endif
endif

# Shared-memory parallelism. libstdc++ routes std::execution::par through TBB,
# which conda-forge does not put in this environment, so a build relying on it
# would quietly run serially. OpenMP ships with the compiler on most toolchains.
#
# Finding the dependency is not the same as being able to link it. On the macOS
# toolchain the flags it supplies drop the default libraries from the link, and
# the failure surfaces as missing unwinder and locale symbols rather than as a
# missing libomp, so probe an actual parallel region before committing to it.
# nvc++ 22.3 trips Meson's OpenMP version sniff (preprocessor
# delimiters). The flag is -mp; do not call dependency('openmp').
openmp_dep_ok = false
_omp_rpath = ''
_omp_libdir = ''
if cpp_id == 'nvidia' or cpp_id == 'nvidia_hpc' or cpp_id == 'pgi'
  openmp_dep = declare_dependency(compile_args: ['-mp'], link_args: ['-mp'])
  openmp_dep_ok = true
  # Host FMA vs device no-FMA (or the reverse) breaks the offload
  # identity check. IEEE division and no FMA on both sides.
  _args += ['-Kieee', '-Mnofma']
else
  openmp_dep = dependency('openmp', required: false)
  openmp_dep_ok = openmp_dep.found()
endif
openmp_found = false
_omp_rpath_dep = []
if openmp_dep_ok
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
    # Clang links libomp by name from its own library directory, such as
    # lib/x86_64-unknown-linux-gnu of an EESSI LLVM or /usr/lib/llvm-21/lib
    # of apt.llvm.org, which the loader does not search; without a runpath
    # there every library and test fails to start unless LD_LIBRARY_PATH
    # names it. A runpath to a directory the loader searches anyway is
    # harmless.
    if cpp_id == 'clang'
      # Linux clang prints libomp.so. Apple clang links libomp.dylib, and
      # -print-file-name=libomp.so then returns the bare name, so the
      # runpath was never recorded and dyld aborted on @rpath/libomp.dylib.
      _omp_names = ['libomp.so']
      if host_machine.system() == 'darwin'
        _omp_names = ['libomp.dylib']
      endif
      foreach _omp_name : _omp_names
        if _omp_libdir == ''
          _omp_file = run_command(cppc.cmd_array() + ['-print-file-name=' + _omp_name],
                                  check: false).stdout().strip()
          if _omp_file.startswith('/') and import('fs').exists(_omp_file)
            _omp_libdir = import('fs').parent(_omp_file)
          endif
        endif
      endforeach
      # conda's clang resource directory does not hold libomp.dylib. The
      # pixi prefix lib directory does, and that is what -fopenmp loads.
      if _omp_libdir == '' and host_machine.system() == 'darwin'
        _conda_omp = run_command('sh', '-c', 'printf %s "${CONDA_PREFIX:-}"',
                                 check: false).stdout().strip()
        if _conda_omp != '' and import('fs').exists(_conda_omp / 'lib' / 'libomp.dylib')
          _omp_libdir = _conda_omp / 'lib'
        endif
      endif
      if _omp_libdir != ''
        _omp_rpath = '-Wl,-rpath,' + _omp_libdir
        _omp_rpath_dep = [declare_dependency(link_args: [_omp_rpath])]
        _deps += _omp_rpath_dep
      endif
    endif
  else
    message('OpenMP was found but a parallel region does not link; ' +
            'building the per-particle loops single-threaded')
  endif
endif

# Snippet and language-front tests compile against this library. They must
# use the same compilers and OpenMP link line. ``cc`` on PATH is not that
# compiler: a pixi environment's cc is conda's gcc, and its glibc 2.28
# sysroot cannot resolve the host compiler's GLIBC and libomp symbols.
_snippet_cc = meson.get_compiler('c').cmd_array()[0]
_snippet_cxx = cppc.cmd_array()[0]
# AddressSanitizer is preloaded so Python can host libseams. That preload
# makes the compiler driver and cc1 fail LeakSanitizer. The wrappers drop
# it and exec the compiler meson used to build the library.
_asan_wrap = get_option('b_sanitize').contains('address')
_real_cc = _snippet_cc
_real_cxx = _snippet_cxx
if _asan_wrap
  _asan_bin = meson.project_source_root() / 'scripts' / 'asan-path'
  _snippet_cc = _asan_bin / 'cc'
  _snippet_cxx = _asan_bin / 'c++'
endif
_snippet_omp = []
if openmp_found
  if cpp_id == 'nvidia' or cpp_id == 'nvidia_hpc' or cpp_id == 'pgi'
    _snippet_omp += ['-mp']
  elif cpp_id == 'clang-cl'
    _snippet_omp += ['-Xclang', '-fopenmp']
  else
    _snippet_omp += ['-fopenmp']
  endif
  if _omp_rpath != ''
    _snippet_omp += [_omp_rpath]
  endif
endif
# A snippet program links against libseams with the compiler meson used,
# and the linker on PATH may be the environment's, whose sysroot is older
# than the host compiler's C and C++ runtimes. libseams then names symbol
# versions that sysroot cannot show at link time; the loader resolves them
# from the host libraries at run time. The library's own executables are
# linked the same way through the environment's LDFLAGS.
if host_machine.system() == 'linux'
  _snippet_omp += ['-Wl,--allow-shlib-undefined']
endif
_snippet_env = environment()
_snippet_env.set('CC', _snippet_cc)
_snippet_env.set('CXX', _snippet_cxx)
if _asan_wrap
  _snippet_env.set('REAL_CC', _real_cc)
  _snippet_env.set('REAL_CXX', _real_cxx)
endif
if _snippet_omp.length() > 0
  _snippet_env.set('SEAMS_OPENMP_LINK', ' '.join(_snippet_omp))
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
    # nvc++ picks a CUDA toolkit that matches the driver. An HPC SDK
    # that does not ship that toolkit fails unless -gpu=cudaX.Y names
    # a bundled release (22.3 ships 11.6; 23.7 ships 12.2).
    offload_try = [
      ['-mp=gpu'],
      ['-mp=gpu', '-gpu=cc80'],
      ['-mp=gpu', '-gpu=cuda12.2'],
      ['-mp=gpu', '-gpu=cc80,cuda12.2'],
      ['-mp=gpu', '-gpu=cuda11.6'],
      ['-mp=gpu', '-gpu=cc80,cuda11.6'],
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
      if cpp_id == 'nvidia' or cpp_id == 'nvidia_hpc' or cpp_id == 'pgi'
        _args += ['-gpu=nofma', '-gpu=nomanaged']
        _link_args += ['-gpu=nofma', '-gpu=nomanaged']
      endif
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
if not hwy_dep.found()
  # Distro packages publish the same library as libhwy.pc.
  hwy_dep = dependency('libhwy', version: '>=1.3.0', required: false)
endif
if hwy_dep.found()
  _deps += [hwy_dep]
  _args += ['-DSEAMS_HAS_HWY=1']
endif

gpulite_opt = get_option('with_gpulite')
gpulite_dep = dependency('gpulite', required: gpulite_opt.disabled() ? gpulite_opt : false,
                         include_type: 'system',
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

# The wrapped linkcell and minimage register their own tests and example
# checks; both follow this project's with_tests.
_sub_flag = get_option('with_tests') ? 'true' : 'false'
_sub_tests = ['with_tests=' + _sub_flag, 'with_examples=' + _sub_flag]
_linkcell_opts = _sub_tests
if gpulite_dep.found() and not gpulite_opt.disabled()
  _linkcell_opts += ['with_gpulite=enabled']
endif
linkcell_dep = dependency('linkcell', version: '>=0.4.2', required: get_option('with_linkcell'),
                          fallback: ['linkcell', 'linkcell_dep'],
                          default_options: _linkcell_opts)
if linkcell_dep.found()
  # Headers only; libyodaLib links the library (src/meson.build).
  _deps += [linkcell_dep.partial_dependency(compile_args: true, includes: true, sources: true)]
  _args += ['-DSEAMS_HAS_LINKCELL=1']
endif

minimage_dep = dependency('minimage', version: '>=0.1.4', required: get_option('with_minimage'),
                          fallback: ['minimage', 'minimage_dep'],
                          default_options: _sub_tests)
if minimage_dep.found()
  _deps += [minimage_dep]
  _args += ['-DSEAMS_HAS_MINIMAGE=1']
endif

chemfiles_libdir = ''
_chem_vendored = false
chemfiles_opt = get_option('with_chemfiles')
# A disabled option must win over a chemfiles.pc on PKG_CONFIG_PATH.
# Otherwise the native frame counters never run.
if chemfiles_opt.disabled()
  chemfiles_dep = dependency('chemfiles-disabled-by-option', required: false)
else
  fs = import('fs')
  chemfiles_dep = dependency('chemfiles', required: false)
  # conda-forge 0.10.4 writes an mmCIF atom with no residue as
  # label_seq_id "." and omits label_entity_id. The reader dereferences
  # that missing column. Reading the written file throws, and a write
  # followed by a read in the same process segfaults inside
  # mmCIFFormat::read. The probe also counts the steps of a MOL2 molecule
  # with an empty counts line, where 0.10.4 segfaults. A library that
  # fails either is replaced by tag 0.10.4 plus cmake/patches/.
  _chemfiles_system_ok = false
  if chemfiles_dep.found()
    _mmcif_probe = meson.project_source_root() / 'cmake' / 'mmcif_no_entity.mmcif'
    _mmcif_run = cppc.run(
      fs.read(meson.project_source_root() / 'cmake' / 'chemfiles_mmcif_probe.cpp'),
      args: ['-DSEAMS_MMCIF_PROBE="' + _mmcif_probe + '"'],
      dependencies: chemfiles_dep,
      name: 'chemfiles reads an mmCIF atom with no residue id and a truncated MOL2')
    _chemfiles_system_ok = _mmcif_run.compiled() and _mmcif_run.returncode() == 0
    if not _chemfiles_system_ok
      message('system chemfiles fails the mmCIF or MOL2 probe; building patched 0.10.4')
      # dependency() has no way to drop a found library. This name is
      # not a package. It clears the rejected system copy so the fetch
      # below is the one that gets linked.
      chemfiles_dep = dependency('chemfiles-mmcif-rejected', required: false)
    endif
  endif
  if not chemfiles_dep.found()
    chem_prefix = meson.project_build_root() / 'chemfiles-prefix'
    cmake_bin = find_program('cmake', required: chemfiles_opt)
    # The prefix records the patches and the CMake files it was built
    # with. A change to either refetches the source, since the patch step
    # runs only on a fetch.
    _chem_patch_ids = []
    foreach _patch : ['patches/chemfiles-0.10.4-mmcif.patch', 'patches/chemfiles-0.10.4-mol2.patch',
                      'patches/chemfiles-0.10.4-frame-visit.patch',
                      'patches/chemfiles-0.10.4-zlib-apple.sh', 'patches/chemfiles-0.10.4-fmt-libcxx.sh',
                      'CMakeLists.txt', 'ChemfilesDep.cmake']
      _chem_patch_ids += [_patch + ' ' + fs.hash(meson.project_source_root() / 'cmake' / _patch, 'sha256')]
    endforeach
    _chem_patch_text = '\n'.join(_chem_patch_ids)
    _chem_stamp = chem_prefix / 'seams-patches.txt'
    _chem_current = fs.exists(chem_prefix / 'lib' / 'libseams-chemfiles.so') and fs.exists(_chem_stamp)
    if _chem_current
      # fs.read() refuses a path in the build tree.
      _chem_current = run_command('cat', _chem_stamp, check: false).stdout().strip() == _chem_patch_text
    endif
    if cmake_bin.found() and not _chem_current
      run_command('rm', '-rf', meson.project_build_root() / 'chemfiles-fetch', chem_prefix,
                  check: true)
      run_command(cmake_bin, '-S', meson.project_source_root() / 'cmake',
                  '-B', meson.project_build_root() / 'chemfiles-fetch',
                  '-DCMAKE_BUILD_TYPE=Release',
                  '-DCMAKE_INSTALL_PREFIX=' + chem_prefix,
                  '-DCMAKE_POLICY_VERSION_MINIMUM=3.5',
                  '-DSEAMS_CHEMFILES_FORCE_VENDORED=ON',
                  check: true)
      run_command(cmake_bin, '--build',
                  meson.project_build_root() / 'chemfiles-fetch',
                  '--parallel', check: true)
      run_command(cmake_bin, '--install',
                  meson.project_build_root() / 'chemfiles-fetch',
                  check: true)
      run_command('sh', '-c', 'printf "%s\\n" "$1" > "$2"', 'sh', _chem_patch_text, _chem_stamp,
                  check: true)
    endif
    # cmake/ChemfilesDep.cmake names the patched library seams-chemfiles,
    # so its soname cannot be satisfied by the unpatched libchemfiles.so.0.10
    # in the conda lib directory the conda driver puts first on the RPATH.
    _chem_real = chem_prefix / 'lib' / 'libseams-chemfiles.so.0.10.4'
    if fs.exists(chem_prefix / 'lib' / 'libseams-chemfiles.so')
      chemfiles_dep = cppc.find_library('seams-chemfiles',
                                        dirs: chem_prefix / 'lib',
                                        required: chemfiles_opt)
      if chemfiles_dep.found()
        chemfiles_dep = declare_dependency(
          compile_args: ['-I' + chem_prefix / 'include', '-DSEAMS_CHEMFILES_PATCHED=1'],
          dependencies: chemfiles_dep)
        chemfiles_libdir = chem_prefix / 'lib'
        _chem_vendored = true
        # The installed libyodaLib finds this copy through its $ORIGIN
        # rpath, under the soname it was linked against.
        if fs.exists(_chem_real)
          install_data(_chem_real, rename: 'libseams-chemfiles.so.0.10',
                       install_dir: get_option('libdir'),
                       install_mode: 'rwxr-xr-x')
        endif
      endif
    endif
  endif
endif
if chemfiles_dep.found()
  _deps += [chemfiles_dep]
  _args += ['-DSEAMS_HAS_CHEMFILES=1']
  probed_libdir = chemfiles_dep.get_variable(pkgconfig: 'libdir',
                                             default_value: '')
  if probed_libdir != ''
    chemfiles_libdir = probed_libdir
  endif
  # The vendored copy sits in the build tree. Meson puts that directory on
  # the build rpath and strips it at install, so only a system chemfiles
  # needs its directory recorded here.
  if chemfiles_libdir != '' and not _chem_vendored
    _link_args += ['-Wl,-rpath,' + chemfiles_libdir]
  endif
endif

# The Rust reader follows the engine's library kind, so a static engine
# embeds it and a shared engine ships it alongside; a shared reader under a
# static engine would point consumers at a build-directory path.
_readcon_kind = get_option('default_library') == 'static' ? 'static' : 'shared'
readcon_opt = get_option('with_readcon')
readcon_dep = dependency('readcon-core', version: '>=0.16.2', required: readcon_opt,
                         fallback: ['readcon-core', 'readcon_dep'],
                         include_type: 'system',
                         default_options: ['default_library=' + _readcon_kind])
if readcon_dep.found()
  _deps += [readcon_dep]
  _args += ['-DSEAMS_HAS_READCON=1']
endif

# Campaign corpora (single LMDB env, or a shards.json root) decode through
# the same readcon-core pin. The db crate vendors its own core for its
# internal parse; seams never frees that handle with this build's
# libreadcon_core.
readcondb_dep = dependency('readcon-db', version: '>=0.1.8',
                           required: readcon_opt.disabled() ? readcon_opt : false,
                           fallback: ['readcon-db', 'readcon_dep'],
                           include_type: 'system',
                           default_options: [
                             'default_library=' + _readcon_kind,
                             'with_tests=false',
                           ])
# The db cdylib sets a soname, so the loader looks up libreadcon_db.so by
# name. The build tree keeps that file in the subproject directory.
_readcondb_rpath = ''
if readcon_dep.found() and readcondb_dep.found()
  _deps += [readcondb_dep]
  _args += ['-DSEAMS_HAS_READCON_DB=1']
  _readcondb_rpath = meson.current_build_dir() / 'subprojects' / 'readcon-db'
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
  ira_inc = meson.project_source_root() / 'src/include/internal'
  if ira_dep.found() and cppc.links(ira_probe, dependencies: [ira_dep],
                                    args: ['-I' + ira_inc],
                                    name: 'libira C API')
    if ira_dep.type_name() == 'pkgconfig'
      ira_libdir = ira_dep.get_variable('libdir', default_value: '')
      if ira_libdir != ''
        ira_dep = declare_dependency(dependencies: [ira_dep],
                                     link_args: ['-Wl,-rpath,' + ira_libdir])
      endif
    endif
    _deps += [ira_dep]
    _args += ['-DSEAMS_HAS_IRA=1']
    ira_found = true
  elif ira_opt.enabled()
    error('with_ira=enabled but libira does not link')
  endif
endif

# spglib (Togo, Shinohara, Tanaka). Optional space groups. The probe
# must run spg_get_dataset, not only find the shared library.
spglib_opt = get_option('with_spglib')
spglib_found = false
if not spglib_opt.disabled()
  spglib_dep = dependency('spglib', required: false)
  if not spglib_dep.found()
    spglib_lib = cppc.find_library('symspg', required: false)
    if spglib_lib.found()
      spglib_dep = declare_dependency(dependencies: [spglib_lib])
    endif
  endif
  spglib_probe = '''
    #include <spglib.h>
    int main() {
      double lattice[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
      double position[1][3] = {{0, 0, 0}};
      int types[1] = {1};
      SpglibDataset *data = spg_get_dataset(lattice, position, types, 1, 1e-4);
      int number = data ? data->spacegroup_number : 0;
      if (data) spg_free_dataset(data);
      return number > 0 ? 0 : 1;
    }
  '''
  if spglib_dep.found() and cppc.links(spglib_probe, dependencies: [spglib_dep],
                                       name: 'spglib dataset')
    if spglib_dep.type_name() == 'pkgconfig'
      spglib_libdir = spglib_dep.get_variable('libdir', default_value: '')
      if spglib_libdir != ''
        spglib_dep = declare_dependency(dependencies: [spglib_dep],
                                        link_args: ['-Wl,-rpath,' + spglib_libdir])
      endif
    endif
    _deps += [spglib_dep]
    _args += ['-DSEAMS_HAS_SPGLIB=1']
    spglib_found = true
  elif spglib_opt.enabled()
    error('with_spglib=enabled but libsymspg does not link')
  endif
endif

# Continuous shape measure (Pinsky, Avnir, Alvarez), GPL-2.0. The
# program is an optional external dependency, off by default, and is
# never vendored into libyodaLib. When enabled, an existing csm_shape
# or shape on PATH is used. Otherwise the program is built from
# third_party/csm if OpenBabel is present. The in-tree measure is the
# published formula on the fixed shells. IRA stays optional under the
# Apache-2.0 option of its dual licence.
csm_shape_bin = ''
csm_opt = get_option('with_csm')
if not csm_opt.disabled()
  csm_prog = find_program('csm_shape', 'shape', required: false)
  if csm_prog.found()
    csm_shape_bin = csm_prog.full_path()
    message('External CSM shape program found (not linked into libyodaLib)')
  else
    openbabel_dep = dependency('openbabel-3', required: false)
    if not openbabel_dep.found()
      openbabel_dep = dependency('openbabel-2.0', required: false)
    endif
    # A module pkg-config can miss the system OpenBabel. The shape program
    # is built by the system compiler against that OpenBabel, and is not
    # linked into libyodaLib.
    if not openbabel_dep.found()
      system_pkg = find_program('/usr/bin/pkg-config', required: false)
      if system_pkg.found()
        foreach ob_name : ['openbabel-3', 'openbabel-2.0']
          ob_probe = run_command(system_pkg, '--exists', ob_name, check: false)
          if ob_probe.returncode() == 0
            openbabel_dep = declare_dependency()
            break
          endif
        endforeach
      endif
    endif
    if not openbabel_dep.found()
      if csm_opt.enabled()
        error('with_csm=enabled but OpenBabel was not found and csm_shape is not on PATH')
      endif
      message('CSM program left off (OpenBabel absent)')
    else
      csm_shape = custom_target(
        'csm_shape',
        input: files(
          'third_party/csm/build.sh',
          'third_party/csm/src/shape.cpp',
          'third_party/csm/src/babelAdapter.cpp',
          'third_party/csm/src/matrix.c',
          'third_party/csm/src/Molecule.c',
          'third_party/csm/src/svdcmp.c',
          'third_party/csm/src/parseFunctions.c',
          'third_party/csm/src/permuter.c',
        ),
        output: 'csm_shape',
        command: [find_program('sh'), files('third_party/csm/build.sh'), '@OUTPUT@', meson.current_source_dir()],
        build_by_default: true,
      )
      csm_shape_bin = csm_shape.full_path()
      message('CSM shape program will be built (not linked into libyodaLib)')
    endif
  endif
endif

sphericart_opt = get_option('with_sphericart')
sphericart_found = false
# Link arguments a program that links libyodaLib, but not sphericart
# itself, needs to resolve libyodaLib's sphericart dependency.
sphericart_consumer_args = []
if not sphericart_opt.disabled()
  sph_lib_dirs = []
  sph_inc_dirs = []
  # The library inside a pip sphericart wheel (often a user site) is only
  # taken when with_sphericart=enabled. In auto mode it made libyodaLib
  # depend on a library that consumers building seams-core as a subproject
  # could not resolve, so auto accepts only pkg-config and CMake packages.
  python_for_sph = find_program('python3', required: false)
  if python_for_sph.found() and sphericart_opt.enabled()
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
        # A cross-prefixed ld (conda's) does not follow libyodaLib's
        # RPATH to libsphericart.so, so a consumer names the folder.
        # -rpath-link is ELF-only. Apple ld records a runtime rpath.
        if host_machine.system() == 'darwin'
          sphericart_consumer_args += ['-Wl,-rpath,' + d]
        else
          sphericart_consumer_args += ['-Wl,-rpath-link,' + d]
        endif
      endforeach
      sph_cflags = []
      foreach d : sph_inc_dirs
        sph_cflags += ['-I' + d]
      endforeach
      sphericart_dep = declare_dependency(dependencies: [sph_lib],
                                          compile_args: sph_cflags,
                                          link_args: sph_rpath + sphericart_consumer_args)
    endif
  endif
  # The probe calls what src/sphericart_ylm.cpp calls, and links without
  # --allow-shlib-undefined (conda's default LDFLAGS carry it), so a
  # library that leaves those symbols or its own dependencies unresolved
  # fails here instead of in a consumer of libyodaLib.
  sphericart_probe = '''
    #include <sphericart.hpp>
    #include <vector>
    int main() {
      sphericart::SphericalHarmonics<double> s(6);
      std::vector<double> xyz{1.0, 0.0, 0.0};
      std::vector<double> out;
      s.compute(xyz, out);
      return out.empty() ? 1 : 0;
    }
  '''
  sph_probe_args = []
  foreach d : sph_inc_dirs
    sph_probe_args += ['-I' + d]
  endforeach
  if cppc.has_link_argument('-Wl,--no-allow-shlib-undefined')
    sph_probe_args += ['-Wl,--no-allow-shlib-undefined']
  endif
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
  'Cutoff fallback (vesin)': vesin_dep.found(),
  'Cutoff list (linkcell pairs_within)': linkcell_dep.found(),
  'Minimum-image wrap (minimage)': minimage_dep.found(),
  'PDB/GRO/DCD readers (chemfiles)': chemfiles_dep.found(),
  '.con reader (readcon-core)': readcon_dep.found(),
  'CON corpus (readcon-db)': readcon_dep.found() and readcondb_dep.found(),
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
  # Locally patched Argum (see third_party/argum/NOTICE.md): a hidden
  # option alias and the -Wconversion/-Wsign-conversion/-Wshadow fixes
  # this project needs are not upstream, so this stays out of Meson
  # WrapDB's dependency() path and out of src/include/internal's
  # clang-tidy and Doxygen scope, like third_party/csm.
  'third_party/argum',
])

# Library
subdir('src')
_snippet_env.set('SEAMS_LIB', seams_api_lib.full_path())

# LAMMPS and GROMACS drivers. They link the C API.
subdir('engines')

nucleation_link = _link_args
if host_system == 'linux'
  nucleation_link += ['-Wl,--enable-new-dtags']
endif

# Installed programs find libyodaLib in the installed libdir. Mach-O
# dylibs carry their absolute install name instead.
_bin_install_rpath = ''
if host_system not in ['darwin', 'windows']
  _bin_install_rpath = '$ORIGIN' / import('fs').relative_to(
    get_option('prefix') / get_option('libdir'),
    get_option('prefix') / get_option('bindir'))
endif
science_front_exe = executable(
  'science-front',
  files('scripts/science_front.cpp') + backward_src,
  dependencies: [yds_dep, backward_dep],
  include_directories: _incdirs,
  cpp_args: _args,
  link_args: nucleation_link,
  build_rpath: meson.current_build_dir() / 'src',
  install_rpath: _bin_install_rpath,
  install: true,
  install_tag: 'bin')

nucleation_pathway_exe = executable(
  'nucleation-pathway',
  files('scripts/nucleation_pathway.cpp') + backward_src,
  dependencies: [yds_dep, backward_dep],
  include_directories: _incdirs,
  cpp_args: _args,
  link_args: nucleation_link,
  build_rpath: meson.current_build_dir() / 'src',
  install_rpath: _bin_install_rpath,
  install: true,
  install_tag: 'bin')

# Committor driver on a 512-molecule mW cell. It skips when LMP is unset.
if get_option('with_tests')
  _pytest = find_program('pytest', required: false)
  if _pytest.found()
    committor_env = environment()
    committor_env.set('NUCLEATION_PATHWAY', nucleation_pathway_exe.full_path())
    committor_env.set('PYTHONDONTWRITEBYTECODE', '1')
    test('committor_driver', _pytest,
      args: ['-q', '-rs', '-p', 'no:cacheprovider',
             meson.project_source_root() / 'tests' / 'test_committor_driver.py'],
      env: committor_env,
      depends: nucleation_pathway_exe,
      suite: 'committor',
      timeout: 1800)
  endif
endif

# Timing driver (bench/analyse.cpp). Outputs stay out of the tree.
subdir('bench')

# Tests
if get_option('with_tests')
  subdir('gpu_host_kernels')
  subdir('tests')
endif

if get_option('with_cli')
  seams_cli_link = _link_args
  if host_machine.system() == 'linux'
    seams_cli_link += ['-Wl,--enable-new-dtags']
  endif
  seams_exe = executable(
    'seams',
    files('src/seams_cli.cpp') + backward_src,
    dependencies: [yds_dep, backward_dep],
    include_directories: _incdirs,
    cpp_args: _args + ['-DSEAMS_VERSION="' + meson.project_version() + '"'],
    link_args: seams_cli_link,
    build_rpath: _readcondb_rpath != '' ? _readcondb_rpath : meson.current_build_dir(),
    install_rpath: _bin_install_rpath,
    install: true)
  brine_front_exe = executable(
    'brine_front',
    files('src/brine_front_main.cpp') + backward_src,
    dependencies: [yds_dep, backward_dep],
    include_directories: _incdirs,
    cpp_args: _args,
    link_args: seams_cli_link,
    build_rpath: _readcondb_rpath != '' ? _readcondb_rpath : meson.current_build_dir(),
    install_rpath: _bin_install_rpath,
    install: true)
endif

# Tests on the CLI. A consumer that builds seams-core as a subproject
# with -Dwith_tests=false registers none of them.
if get_option('with_cli') and get_option('with_tests')
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
  # In-situ engine runs checked frame by frame against this CLI.
  subdir('engines/e2e')

  python3 = find_program('python3', required: true)
  test(
    'docs_language_examples',
    python3,
    args: [meson.project_source_root() / 'docs' / 'check_language_examples.py'],
    suite: 'docs',
    timeout: 60)
  test(
    'docs_c_api_tour',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'check_c_api_tour.py',
      meson.project_source_root(),
      meson.project_build_root() / 'src',
      meson.project_build_root(),
    ],
    env: _snippet_env,
    suite: 'docs',
    timeout: 120)
  test(
    'docs_check',
    python3,
    args: [meson.project_source_root() / 'docs' / 'check_examples.py'],
    suite: 'docs',
    timeout: 60)
  test(
    'docs_table_csv',
    python3,
    args: [meson.project_source_root() / 'docs' / 'check_tables.py'],
    suite: 'docs',
    timeout: 60)
  test(
    'docs_prose',
    python3,
    args: [meson.project_source_root() / 'docs' / 'check_prose.py'],
    suite: 'docs',
    timeout: 60)
  test(
    'release_versions',
    python3,
    args: [meson.project_source_root() / 'scripts' / 'check_release_version.py'],
    suite: 'release',
    timeout: 60)
  test(
    'release_versions_selftest',
    python3,
    args: ['-m', 'doctest', meson.project_source_root() / 'scripts' / 'check_release_version.py'],
    suite: 'release',
    timeout: 60)
  test(
    'docs_llms_txt',
    python3,
    args: [meson.project_source_root() / 'scripts' / 'gen_llms_txt.py', '--check'],
    suite: 'docs',
    timeout: 60)
  test(
    'docs_llms_txt_selftest',
    python3,
    args: ['-m', 'doctest', meson.project_source_root() / 'scripts' / 'gen_llms_txt.py'],
    workdir: meson.project_source_root() / 'scripts',
    suite: 'docs',
    timeout: 60)
  test(
    'docs_discoverability',
    python3,
    args: [meson.project_source_root() / 'docs' / 'check_discoverability.py', '--self-test'],
    suite: 'docs',
    timeout: 60)
  test(
    'docs_c_symbols',
    python3,
    args: ['-m', 'doctest', meson.project_source_root() / 'docs' / 'check_c_symbols.py'],
    workdir: meson.project_source_root() / 'docs',
    suite: 'docs',
    timeout: 60)
  test(
    'docs_api_coverage',
    python3,
    args: [meson.project_source_root() / 'docs' / 'check_api_coverage.py'],
    suite: 'docs',
    timeout: 60)
  test(
    'docs_redirects',
    python3,
    args: [meson.project_source_root() / 'docs' / 'check_redirects.py'],
    suite: 'docs',
    timeout: 60)
  test(
    'docs_binding_api',
    python3,
    args: [meson.project_source_root() / 'docs' / 'check_binding_api.py'],
    suite: 'docs',
    timeout: 300)
  tutorial_env = environment()
  tutorial_env.prepend('LD_LIBRARY_PATH', meson.project_build_root() / 'src')
  tutorial_env.set('SEAMS_LIB', seams_api_lib.full_path())
  # bindings/ is not configured yet. The filename is the one shared_library()
  # will write; the loader uses it only when that file exists.
  _lang_file = 'libseams_lang.so'
  if host_machine.system() == 'darwin'
    _lang_file = 'libseams_lang.dylib'
  elif host_machine.system() == 'windows'
    _lang_file = 'libseams_lang.dll'
  endif
  tutorial_env.set('SEAMS_LANG_LIB', meson.project_build_root() / 'bindings' / _lang_file)
  tutorial_env.set('CC', _snippet_cc)
  tutorial_env.set('CXX', _snippet_cxx)
  if _asan_wrap
    tutorial_env.set('REAL_CC', _real_cc)
    tutorial_env.set('REAL_CXX', _real_cxx)
  endif
  if _snippet_omp.length() > 0
    tutorial_env.set('SEAMS_OPENMP_LINK', ' '.join(_snippet_omp))
  endif
  test(
    'docs_tutorial_clathrate',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'tutorials' / 'test_clathrate.py',
      meson.project_source_root(),
      seams_exe.full_path(),
      meson.project_build_root() / 'src',
    ],
    depends: [seams_exe],
    suite: 'docs',
    timeout: 180)
  test(
    'docs_tutorial_interface',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'tutorials' / 'test_interface.py',
      meson.project_source_root(),
      seams_exe.full_path(),
      meson.project_build_root() / 'src',
    ],
    env: tutorial_env,
    depends: [seams_exe],
    suite: 'docs',
    timeout: 120)
  foreach _paper : ['interface', 'cages', 'rings', 'brine', 'nucleation']
    test(
      'docs_paper_' + _paper,
      python3,
      args: [
        meson.project_source_root() / 'docs' / 'tutorials' / 'test_paper.py',
        _paper,
        meson.project_source_root(),
        seams_exe.full_path(),
        meson.project_build_root() / 'src',
        meson.project_build_root(),
      ],
      env: tutorial_env,
      depends: [seams_exe],
      suite: 'docs',
      timeout: 300)
  endforeach
  # Full paper runs. Suite slow only, so `meson test --suite docs` stays fast.
  # Off unless -Dpaper_full=true. The script then exits 77 and says so.
  paper_env = environment()
  paper_env.prepend('LD_LIBRARY_PATH', meson.project_build_root() / 'src')
  if get_option('paper_full')
    paper_env.set('SEAMS_PAPER_FULL', '1')
  endif
  foreach _paper : ['interface', 'cages', 'rings', 'brine', 'nucleation']
    test(
      'docs_paper_' + _paper + '_full',
      python3,
      args: [
        meson.project_source_root() / 'docs' / 'tutorials' / 'test_paper_full.py',
        _paper,
        meson.project_source_root(),
        seams_exe.full_path(),
        meson.project_build_root() / 'src',
        meson.project_build_root(),
      ],
      env: paper_env,
      depends: [seams_exe],
      suite: 'slow',
      timeout: 900)
  endforeach
  test(
    'docs_frame_tutorials',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'tutorials' / 'test_hf_frames.py',
      '--small',
      meson.project_source_root(),
      seams_exe.full_path(),
      meson.project_build_root() / 'src',
    ],
    env: tutorial_env,
    depends: [seams_exe, seams_api_lib],
    suite: 'docs',
    timeout: 180)
  test(
    'docs_tutorial_nucleation',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'tutorials' / 'test_nucleation.py',
      meson.project_source_root(),
      meson.project_build_root() / 'src',
    ],
    env: tutorial_env,
    depends: [seams_api_lib],
    suite: 'docs',
    timeout: 180)
  test(
    'docs_tutorial_brine',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'tutorials' / 'test_brine.py',
      meson.project_source_root(),
      meson.project_build_root() / 'src',
    ],
    env: tutorial_env,
    depends: [seams_api_lib],
    suite: 'docs',
    timeout: 120)
  test(
    'docs_science_front',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'howto' / 'test_science_front.py',
      meson.project_source_root(),
      meson.project_build_root(),
    ],
    env: tutorial_env,
    suite: 'docs',
    timeout: 180)
  test(
    'docs_tutorial_insitu',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'tutorials' / 'test_insitu.py',
      meson.project_source_root(),
      meson.project_build_root(),
    ],
    env: tutorial_env,
    suite: 'docs',
    timeout: 180)
  test(
    'engines_tutorial',
    python3,
    args: [
      meson.project_source_root() / 'engines' / 'tutorial_walk.py',
      meson.project_build_root(),
    ],
    env: tutorial_env,
    depends: [seams_api_lib],
    suite: 'docs',
    timeout: 120)
  test(
    'seams_water_types_cutoff',
    python3,
    args: [
      meson.project_source_root() / 'tests' / 'test_water_mask.py',
      seams_exe.full_path(),
    ],
    depends: [seams_exe],
    timeout: 60)
  test(
    'seams_mol_rings',
    python3,
    args: [
      meson.project_source_root() / 'tests' / 'test_mol_rings.py',
      seams_exe.full_path(),
    ],
    depends: [seams_exe],
    timeout: 60)

  test(
    'check_option_names',
    python3,
    args: [
      meson.project_source_root() / 'tests' / 'check_option_names.py',
      meson.project_source_root(),
      seams_exe.full_path(),
    ],
    depends: [seams_exe],
    timeout: 60)
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
    'seams_chill_plus_default_is_cutoff_list',
    python3,
    args: [
      '-c',
      'import subprocess,sys; run=lambda *a: subprocess.run([sys.argv[1],"chill-plus",sys.argv[2],"--cutoff","3.5",*a],capture_output=True,text=True); d=run(); m=run("--graph","mutual"); ok=d.returncode==0 and m.returncode==0 and "interClathrate 9 water 241" in d.stdout and "clathrate 1 interClathrate 7 water 242" in m.stdout; sys.exit(0 if ok else 1)',
      seams_exe.full_path(),
      meson.project_source_root() / 'input' / 'traj' / 'exampleTraj.lammpstrj',
    ],
    depends: [seams_exe])
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
      'import subprocess,sys; p=subprocess.run([sys.argv[1],"cn",sys.argv[2],"--pair-types","2,2","--cutoff","3.5"],capture_output=True,text=True); text=p.stdout+p.stderr; sys.exit(0 if p.returncode==0 and "site-site" in text else 1)',
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
  test(
    'seams_carried_frames_match_rebuild',
    python3,
    args: [
      files('tests/cli_carry.py'),
      seams_exe.full_path(),
      meson.project_source_root() / 'input' / 'traj',
    ],
    depends: [seams_exe],
    timeout: 600)
  _mpiexec = find_program('mpiexec', required: false)
  if mpi_found and _mpiexec.found()
    test(
      'seams_mpi_frames_match_serial',
      python3,
      args: [
        files('tests/mpi_frames.py'),
        _mpiexec.full_path(),
        seams_exe.full_path(),
        meson.project_source_root() / 'input' / 'traj' / 'mW_cubic.lammpstrj',
        meson.project_source_root() / 'input' / 'eon' / 'tiny_multi_cuh2.con',
      ],
      depends: [seams_exe],
      is_parallel: false,
      timeout: 300)
  endif
endif

subdir('bindings')
_snippet_env.set('SEAMS_LANG_LIB', seams_lang.full_path())

if get_option('with_cli') and get_option('with_tests')
  # check_installed.py installs with --no-rebuild, because a rebuild here
  # would relink libyodaLib and libseams under the tests that load them.
  # Every installed target is therefore built before the test starts.
  _installed_targets = [seams_api_lib, ydslib, seams_lang, seams_exe, brine_front_exe,
                        science_front_exe, nucleation_pathway_exe] + backward_cpp_installed
  foreach _optional : ['run_hook', 'lammps_plugin', 'gmx_seams', 'openmm_plugin', 'seams_plumed']
    if is_variable(_optional)
      _installed_targets += get_variable(_optional)
    endif
  endforeach
  test(
    'c_api_installed',
    python3,
    args: [
      meson.project_source_root() / 'examples' / 'c' / 'check_installed.py',
      meson.project_build_root(),
      meson.project_source_root(),
    ],
    env: _snippet_env,
    suite: 'docs',
    depends: _installed_targets,
    timeout: 600)
  # ldd and LD_BIND_NOW are the ELF loader's.
  if host_machine.system() == 'linux'
    test(
      'seams_lang_pkgconfig',
      python3,
      args: [
        meson.project_source_root() / 'tests' / 'check_lang_pkgconfig.py',
        meson.project_build_root(),
        meson.project_source_root(),
      ],
      env: _snippet_env,
      depends: _installed_targets,
      timeout: 600)
  endif
  test(
    'tests_do_not_rebuild',
    python3,
    args: [
      meson.project_source_root() / 'tests' / 'check_no_rebuild.py',
      meson.project_source_root(),
      meson.project_build_root(),
    ],
    timeout: 120)
  test(
    'docs_migrate',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'tutorials' / 'test_migrate.py',
      meson.project_source_root(),
      seams_exe.full_path(),
      meson.project_build_root() / 'src',
      meson.project_build_root() / 'bindings',
    ],
    env: tutorial_env,
    depends: [seams_exe, seams_api_lib, seams_lang],
    suite: 'docs',
    timeout: 180)
  test(
    'docs_chill_rule',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'test_chill_rule_page.py',
      meson.project_source_root(),
      seams_exe.full_path(),
      meson.project_build_root() / 'src',
      meson.project_build_root() / 'bindings',
      meson.project_build_root(),
    ],
    env: tutorial_env,
    depends: [seams_exe, seams_api_lib, seams_lang],
    suite: 'docs',
    timeout: 120)
  molecule_env = environment()
  molecule_env.prepend('LD_LIBRARY_PATH', meson.project_build_root() / 'src')
  # Full publication stills are 4096 px and 256 samples. The test checks
  # the frame and the palette, which a 2048 px, 64-sample draw already shows.
  molecule_env.set('MOLECULES_SCALE', '1')
  molecule_env.set('MOLECULES_SAMPLES', '64')
  test(
    'docs_render_molecules',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'test_render_molecules.py',
      '--slug', 'ice-ih',
      '--out', meson.project_build_root() / 'molecule-stills',
      '--libdir', meson.project_build_root() / 'src',
    ],
    env: molecule_env,
    depends: [seams_api_lib],
    suite: 'docs',
    timeout: 900)
  test(
    'docs_video_cards',
    python3,
    args: [meson.project_source_root() / 'docs' / 'check_video_cards.py'],
    suite: 'docs',
    timeout: 60)
  test(
    'docs_topic_embeds',
    python3,
    args: [meson.project_source_root() / 'docs' / 'check_topic_embeds.py'],
    suite: 'docs',
    timeout: 60)
  test(
    'docs_video_tapes',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'video_tapes.py',
      meson.project_source_root(),
      seams_exe.full_path(),
      meson.project_build_root() / 'src',
      meson.project_build_root() / 'bindings',
      meson.project_build_root(),
    ],
    env: tutorial_env,
    depends: [seams_exe, seams_api_lib, seams_lang],
    suite: 'docs',
    timeout: 600)
  test(
    'docs_gallery',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'render_gallery.py',
      '--check',
      seams_exe.full_path(),
      meson.project_build_root() / 'src',
      meson.project_build_root() / 'bindings',
    ],
    env: tutorial_env,
    depends: [seams_exe, seams_api_lib, seams_lang],
    suite: 'docs',
    timeout: 300)
  binding_env = environment()
  binding_env.prepend('LD_LIBRARY_PATH', meson.project_build_root() / 'src')
  binding_env.prepend('LD_LIBRARY_PATH', meson.project_build_root() / 'bindings')
  binding_env.set('SEAMS_LIB', seams_api_lib.full_path())
  binding_env.set('SEAMS_LANG_LIB', seams_lang.full_path())
  binding_env.set('CC', _snippet_cc)
  binding_env.set('CXX', _snippet_cxx)
  if _asan_wrap
    binding_env.set('REAL_CC', _real_cc)
    binding_env.set('REAL_CXX', _real_cxx)
  endif
  if _snippet_omp.length() > 0
    binding_env.set('SEAMS_OPENMP_LINK', ' '.join(_snippet_omp))
  endif
  test(
    'docs_lang_ice',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'examples' / 'run_lang_ice.py',
      meson.project_source_root(),
      meson.project_build_root(),
    ],
    env: binding_env,
    depends: [seams_api_lib, seams_lang],
    suite: 'docs',
    timeout: 120)
  test(
    'docs_science_render_copy',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'check_science_render.py',
    ],
    suite: 'docs',
    timeout: 60)
  test(
    'docs_interface_layers',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'check_interface_layers.py',
    ],
    env: binding_env,
    depends: [seams_api_lib],
    suite: 'docs',
    timeout: 60)
  test(
    'docs_figure_layout',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'check_figure_layout.py',
    ],
    suite: 'docs',
    timeout: 60)
  test(
    'docs_brine_cli',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'check_brine_cli.py',
    ],
    suite: 'docs',
    timeout: 180)
  test(
    'docs_figure_census',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'check_figure_census.py',
    ],
    suite: 'docs',
    timeout: 60)
  test(
    'docs_palette',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'check_palette.py',
    ],
    suite: 'docs',
    timeout: 60)
  test(
    'docs_cage_edges',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'check_cage_edges.py',
    ],
    env: binding_env,
    depends: [seams_api_lib],
    suite: 'docs',
    timeout: 120)
  test(
    'docs_engine_pages',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'check_engine_pages.py',
    ],
    suite: 'docs',
    timeout: 60)
  test(
    'docs_sh_cell',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'check_sh_cell.py',
    ],
    env: binding_env,
    depends: [seams_lang],
    suite: 'docs',
    timeout: 60)
  test(
    'docs_science_tutorials',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'examples' / 'science_tutorials.py',
      meson.project_source_root(),
      meson.project_build_root(),
    ],
    env: binding_env,
    depends: [seams_api_lib],
    suite: 'docs',
    # ASan makes the render CPU-bound; keep it from sharing the runner cores.
    is_parallel: not get_option('b_sanitize').contains('address'),
    timeout: 180)
  test(
    'docs_engine_walkthrough',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'examples' / 'engine_walkthrough.py',
      meson.project_source_root(),
      meson.project_build_root(),
    ],
    env: binding_env,
    depends: [seams_api_lib],
    suite: 'docs',
    timeout: 120)
  test(
    'docs_capi_symbols',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'check_capi_symbols.py',
    ],
    suite: 'docs',
    timeout: 60)
  test(
    'docs_capi_examples',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'examples' / 'run_capi_examples.py',
      meson.project_source_root(),
      meson.project_build_root(),
    ],
    env: binding_env,
    depends: [seams_api_lib],
    suite: 'docs',
    timeout: 120)
  test(
    'docs_core_counts',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'examples' / 'run_core_counts.py',
      meson.project_source_root(),
      meson.project_build_root(),
    ],
    env: binding_env,
    depends: [seams_api_lib, seams_lang],
    suite: 'docs',
    timeout: 300)
  test(
    'docs_examples',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'examples' / 'run_examples.py',
      meson.project_source_root(),
      meson.project_build_root(),
    ],
    env: binding_env,
    depends: [seams_api_lib, seams_lang, seams_exe],
    suite: 'docs',
    timeout: 300)
  foreach _engine : ['gromacs', 'lammps', 'openmm', 'plumed', 'insitu']
    test(
      'docs_engine_' + _engine,
      python3,
      args: [
        meson.project_source_root() / 'docs' / 'howto' / 'test_engines.py',
        _engine,
        meson.project_source_root(),
        meson.project_build_root(),
      ],
      env: binding_env,
      depends: [seams_lang],
      suite: 'docs',
      timeout: 180)
  endforeach
  test(
    'docs_doc_blocks',
    python3,
    args: [
      meson.project_source_root() / 'docs' / 'run_doc_blocks.py',
      meson.project_build_root(),
    ],
    env: binding_env,
    depends: [seams_exe, seams_api_lib, seams_lang],
    suite: 'docs',
    timeout: 600)
  foreach _lang : ['c', 'cpp']
    test(
      'docs_lang_' + _lang,
      python3,
      args: [
        meson.project_source_root() / 'docs' / 'examples' / 'run_lang_tasks.py',
        _lang,
        meson.project_build_root(),
      ],
      env: binding_env,
      depends: [seams_api_lib, seams_lang],
      suite: 'docs',
      timeout: 180)
  endforeach
  _snippet_list = run_command(
    python3,
    meson.project_source_root() / 'docs' / 'run_snippets.py',
    '--list',
    check: true).stdout().strip().split('\n')
  foreach _page : _snippet_list
    if _page != ''
      test(
        'docs_snippet_' + _page,
        python3,
        args: [
          meson.project_source_root() / 'docs' / 'run_snippets.py',
          _page,
          meson.project_build_root(),
        ],
        env: binding_env,
        depends: [seams_api_lib, seams_lang, binding_oracle],
        suite: 'docs',
        timeout: 300)
    endif
  endforeach
  test(
    'docs_snippet_harness',
    python3,
    args: ['-m', 'doctest', meson.project_source_root() / 'docs' / 'run_snippets.py'],
    workdir: meson.project_source_root() / 'docs',
    suite: 'docs',
    timeout: 60)
endif

if get_option('with_tests')
  # C API proof. Sanitizers and valgrind are separate configure options.
  # clang-cl rejects a raw -std=c11. override_options lets meson emit /std:c11.
  # libm is in the C runtime on Windows.
  _libm = host_machine.system() == 'windows' ? [] : ['-lm']
  capi_suite = executable(
    'capi_suite',
    files('tests/capi_suite.c'),
    link_with: seams_api_lib,
    include_directories: include_directories('include'),
    override_options: ['c_std=c11'],
    c_args: ['-Wall', '-Wextra', '-Werror', '-I' + meson.current_build_dir() / 'src'],
    link_args: _libm)
  capi_cpp_ref = executable(
    'capi_cpp_ref',
    files('tests/capi_cpp_ref.cpp'),
    dependencies: [yds_dep],
    include_directories: _incdirs,
    cpp_args: _args,
    link_args: _link_args,
    build_rpath: meson.current_build_dir() / 'src')
  _capi_env = environment()
  _capi_env.prepend('LD_LIBRARY_PATH', meson.current_build_dir() / 'src')
  if _capi_san == 'thread'
    # libgomp here is not TSan-instrumented, so the barrier between two
    # parallel loops is invisible and reported as a race. The suite still
    # fails when q6 at 1, 2, 4 and 8 threads disagrees.
    _capi_env.set('TSAN_OPTIONS', 'halt_on_error=0:exitcode=0')
  endif
  if _capi_san.contains('address')
    _capi_env.set('LSAN_OPTIONS',
                  'suppressions=' + meson.project_source_root() / 'fuzz' / 'lsan.supp')
  endif
  _capi_args = [
    meson.project_source_root(),
    capi_cpp_ref.full_path(),
    meson.project_source_root() / 'tests' / 'capi_ovito_freud.py',
    meson.current_build_dir() / 'capi-work',
  ]
  _capi_run = capi_suite
  _capi_run_args = _capi_args
  if _capi_san == 'thread' and host_machine.system() == 'linux'
    # GCC's libtsan expects the 28-bit mmap randomisation of older kernels.
    # Under 32 bits it finds the executable inside its shadow range and
    # stops with "unexpected memory mapping". setarch -R turns address
    # randomisation off for the suite and the capi_cpp_ref it starts.
    _capi_run = find_program('setarch', required: true)
    _capi_run_args = [host_machine.cpu(), '-R', capi_suite.full_path()] + _capi_args
  endif
  test(
    'capi_suite',
    _capi_run,
    args: _capi_run_args,
    env: _capi_env,
    depends: [capi_suite, capi_cpp_ref, seams_api_lib],
    timeout: 3600,
    is_parallel: false)
  _cross_env = environment()
  _cross_env.prepend('LD_LIBRARY_PATH', meson.current_build_dir() / 'src')
  _cross_env.set('SEAMS_LIB', seams_api_lib.full_path())
  _cross_env.set('OMP_NUM_THREADS', '1')
  _cross_env.set('RAYON_NUM_THREADS', '1')
  _cross_env.set('OVITO_THREAD_COUNT', '1')
  test(
    'capi_crosscheck',
    find_program('python3', required: true),
    args: [
      meson.project_source_root() / 'tests' / 'capi_ovito_freud.py',
      '--root', meson.project_source_root(),
    ],
    env: _cross_env,
    depends: [seams_api_lib],
    timeout: 900,
    is_parallel: false)
  if get_option('capi_valgrind')
    test(
      'capi_valgrind',
      find_program('valgrind', required: true),
      args: [
        '--leak-check=full', '--show-leak-kinds=definite',
        '--errors-for-leak-kinds=definite', '--error-exitcode=1',
        capi_suite.full_path(),
      ] + _capi_args,
      env: _capi_env,
      depends: [capi_suite, capi_cpp_ref],
      timeout: 7200,
      is_parallel: false)
  endif

  _capi_warn = ['-Wall', '-Wextra', '-Wpedantic', '-Werror']
  _capi_headers = [
    ['seams_h', 'seams.h',
     ['-I' + meson.project_source_root() / 'include',
      '-I' + meson.current_build_dir() / 'src']],
    ['seams_version_h', 'seams_version.h',
     ['-I' + meson.current_build_dir() / 'src']],
    ['seams_lang_h', 'seams_lang.h',
     ['-I' + meson.project_source_root() / 'bindings' / 'cxx']],
  ]
  foreach _hdr : _capi_headers
    test('capi_header_c99_' + _hdr[0],
      executable('capi_header_c99_' + _hdr[0], files('tests/capi_header_unit.c'),
        override_options: ['c_std=c99'],
        c_args: ['-include', _hdr[1]] + _capi_warn + _hdr[2]))
    test('capi_header_cxx17_' + _hdr[0],
      executable('capi_header_cxx17_' + _hdr[0], files('tests/capi_header_unit.cpp'),
        override_options: ['cpp_std=c++17'],
        cpp_args: ['-include', _hdr[1]] + _capi_warn + _hdr[2]))
  endforeach

  capi_contract = executable(
    'capi_contract',
    files('tests/capi_contract.c'),
    link_with: seams_api_lib,
    include_directories: include_directories('include'),
    override_options: ['c_std=c11'],
    c_args: _capi_warn + ['-I' + meson.current_build_dir() / 'src'],
    link_args: _libm)
  test(
    'capi_contract',
    capi_contract,
    args: [
      meson.project_source_root() / 'tests' / 'capi_symbols.txt',
      meson.project_source_root(),
    ],
    env: _capi_env,
    timeout: 300)
  capi_thread_scale = executable(
    'capi_thread_scale',
    files('tests/capi_thread_scale.c'),
    link_with: ydslib,
    dependencies: [dependency('openmp', language: 'c', required: false)] + _omp_rpath_dep,
    include_directories: include_directories('src/include/internal'),
    c_args: ['-std=c11'] + _capi_warn,
    link_args: ['-lm'],
    build_rpath: meson.current_build_dir() / 'src')
  # Alone, so the logged times are not shared with other tests.
  test(
    'capi_thread_scale',
    capi_thread_scale,
    args: [meson.current_build_dir() / 'capi_thread_scale.txt'],
    env: _capi_env,
    timeout: 1800,
    is_parallel: false)
  foreach _family : ['neighbours', 'rings', 'chill_plus', 'q4q6', 'cages', 'errors']
    test(
      'capi_example_' + _family,
      executable(
        'capi_example_' + _family,
        files('docs/examples/c/' + _family + '.c'),
        link_with: seams_api_lib,
        include_directories: include_directories('include'),
        override_options: ['c_std=c11'],
        c_args: _capi_warn + ['-I' + meson.current_build_dir() / 'src']),
      args: _family == 'errors' ? [] : [meson.project_source_root()],
      env: _capi_env,
      timeout: 300)
  endforeach
  _installed_examples_env = {'CC': ' '.join(meson.get_compiler('c').cmd_array()),
                             'SEAMS_OPENMP_LINK': ' '.join(_snippet_omp)}
  if _asan_wrap
    _installed_examples_env = {'CC': _snippet_cc, 'REAL_CC': _real_cc,
                               'SEAMS_OPENMP_LINK': ' '.join(_snippet_omp)}
  endif
  test(
    'capi_installed_examples',
    find_program('python3', required: true),
    args: [
      meson.project_source_root() / 'docs' / 'examples' / 'c' / 'check_installed.py',
      meson.project_build_root(),
      meson.project_source_root(),
    ],
    env: _installed_examples_env,
    depends: [seams_api_lib, ydslib],
    timeout: 600)
  test(
    'capi_const_threads',
    executable(
      'capi_const_threads',
      files('tests/capi_const_threads.c'),
      link_with: seams_api_lib,
      include_directories: include_directories('include'),
      override_options: ['c_std=c11'],
      c_args: _capi_warn + ['-I' + meson.current_build_dir() / 'src'],
      dependencies: dependency('threads'),
      link_args: _libm),
    env: _capi_env,
    timeout: 300)
  # GCC marks only ASan and TSan builds; clang reports every sanitizer.
  _signal_mode = (_capi_san.contains('address') or _capi_san.contains('thread') or
                  (cppc.get_id() == 'clang' and _capi_san != 'none')) ? 'sanitized' : 'plain'
  test(
    'capi_signal_owner',
    executable(
      'capi_signal_owner',
      files('tests/capi_signal_owner.c'),
      link_with: seams_api_lib,
      include_directories: include_directories('include'),
      override_options: ['c_std=c11'],
      c_args: _capi_warn + ['-I' + meson.current_build_dir() / 'src']),
    args: [_signal_mode],
    env: _capi_env,
    timeout: 60)
  test(
    'capi_q_threads',
    executable(
      'capi_q_threads',
      files('tests/capi_q_threads.c'),
      link_with: seams_api_lib,
      include_directories: include_directories('include'),
      override_options: ['c_std=c11'],
      c_args: _capi_warn + ['-I' + meson.current_build_dir() / 'src'],
      dependencies: dependency('threads'),
      link_args: _libm),
    env: _capi_env,
    timeout: 300)
  executable(
    'bench_capi_rss',
    files('tests/bench_capi_rss.c'),
    link_with: seams_api_lib,
    include_directories: include_directories('include'),
    override_options: ['c_std=c11'],
    c_args: _capi_warn + ['-I' + meson.current_build_dir() / 'src'],
    link_args: _libm)
  test(
    'capi_symbols',
    find_program('python3', required: true),
    args: [
      meson.project_source_root() / 'tests' / 'capi_symbols.py',
      meson.project_source_root() / 'include' / 'seams.h',
      meson.project_source_root() / 'tests' / 'capi_symbols.txt',
      seams_api_lib.full_path(),
    ],
    depends: [seams_api_lib],
    timeout: 60)
  test(
    'capi_symbols_selftest',
    find_program('python3', required: true),
    args: [
      meson.project_source_root() / 'tests' / 'capi_symbols_selftest.py',
      meson.project_source_root() / 'include' / 'seams.h',
      meson.project_source_root() / 'tests' / 'capi_symbols.txt',
      seams_api_lib.full_path(),
    ],
    depends: [seams_api_lib],
    timeout: 60)
  _capi_docs_args = [
    meson.project_source_root() / 'include' / 'seams.h',
    meson.project_source_root() / 'tests' / 'capi_symbols.txt',
    meson.project_source_root() / 'docs' / 'api' / 'c.md',
  ]
  test(
    'capi_docs_complete',
    find_program('python3', required: true),
    args: [meson.project_source_root() / 'docs' / 'capi_docs.py'] + _capi_docs_args,
    timeout: 60)
  test(
    'capi_docs_selftest',
    find_program('python3', required: true),
    args: [meson.project_source_root() / 'docs' / 'capi_docs_selftest.py'] + _capi_docs_args,
    timeout: 60)
  _doxygen = find_program('doxygen', required: false)
  if _doxygen.found()
    _doxy_data = configuration_data()
    _doxy_data.set('VERSION', meson.project_version())
    _doxy_data.set('SRCDIR', meson.project_source_root())
    _doxy_data.set('BUILDDIR', meson.project_build_root())
    _doxy_data.set('OUTDIR', meson.current_build_dir() / 'capi-doxygen')
    _doxyfile = configure_file(
      input: 'docs/Doxyfile.capi.in',
      output: 'Doxyfile.capi',
      configuration: _doxy_data)
    test(
      'capi_doxygen',
      _doxygen,
      args: [_doxyfile],
      workdir: meson.current_build_dir(),
      timeout: 120)
  endif

  _abi_env = environment()
  _abi_depends = [seams_api_lib]
  if _abidiff.found()
    _abi_env.set('SEAMS_ABIDIFF', _abidiff.full_path())
    _abi_env.set('SEAMS_ABI_DEBUG_LIB', seams_abi_lib.full_path())
    _abi_depends += seams_abi_lib
  endif
  test(
    'capi_abi',
    find_program('sh'),
    args: [
      meson.project_source_root() / 'scripts' / 'abi_check.sh',
      meson.current_build_dir(),
      meson.project_source_root() / 'scripts' / 'abi_baseline.txt',
    ],
    env: _abi_env,
    depends: _abi_depends,
    timeout: 300)
  _lang_abi_env = environment()
  _lang_abi_depends = []
  if _abidiff.found()
    _lang_abi_env.set('SEAMS_ABIDIFF', _abidiff.full_path())
    _lang_abi_depends += seams_lang_abi_lib
  endif
  foreach _mode : ['check', 'selftest']
    test(
      _mode == 'check' ? 'capi_abi_2_13' : 'capi_abi_2_13_selftest',
      find_program('python3', required: true),
      args: [
        meson.project_source_root() / 'scripts' / 'abi_check_lang.py',
        _mode,
        meson.project_source_root() / 'scripts' / 'abi_2.13.0_seams_lang.abi',
        meson.project_source_root() / 'scripts' / 'abi_seams_lang.supp',
        meson.project_source_root() / 'bindings' / 'cxx',
        meson.current_build_dir() / 'bindings' / 'abi' / 'libseams_lang.so',
        meson.current_build_dir() / 'abi-2.13-' + _mode,
      ],
      env: _lang_abi_env,
      depends: _lang_abi_depends,
      timeout: 300)
  endforeach

  # libFuzzer harnesses. Under clang they are built with -fsanitize=fuzzer,
  # plus the capi_sanitize sanitizer, so libFuzzer gets coverage from the
  # harness; capi_fuzz_coverage instruments libseams as well. Other
  # compilers link a clang fuzzer runtime from the usual Debian path,
  # without coverage, when one is installed.
  _fuzz_cc = meson.get_compiler('c')
  _fuzz_args = []
  if _fuzz_cc.get_id() == 'clang' and _fuzz_cc.has_argument('-fsanitize=fuzzer-no-link')
    _fuzz_args = ['-fsanitize=fuzzer' + (_capi_san == 'none' ? '' : ',' + _capi_san)]
  elif _capi_san == 'none'
    foreach _cand : [
      '/usr/lib/llvm-18/lib/clang/18/lib/linux/libclang_rt.fuzzer-x86_64.a',
      '/usr/lib/llvm-19/lib/clang/19/lib/linux/libclang_rt.fuzzer-x86_64.a',
      '/usr/lib/llvm-20/lib/clang/20/lib/linux/libclang_rt.fuzzer-x86_64.a',
    ]
      if _fuzz_args.length() == 0 and import('fs').exists(_cand)
        _fuzz_args = [_cand]
      endif
    endforeach
  endif
  if _fuzz_args.length() > 0
    _fuzz_env = environment()
    _fuzz_env.set('LSAN_OPTIONS',
                  'suppressions=' + meson.project_source_root() / 'fuzz' / 'lsan.supp')
    # The compiler invokes whichever `ld` is first on PATH. A conda prefix
    # puts its sysroot ld there. That ld links the prefix's libpthread and
    # libdl, whose GLIBC_PRIVATE symbols are not in the host libc this
    # compiler was built against. -B names the host linker when the
    # compiler itself is not the conda one.
    _fuzz_bindir = []
    _conda_prefix = run_command('sh', '-c', 'printf %s "${CONDA_PREFIX:-}"', check: false).stdout().strip()
    _cc_path = find_program(_fuzz_cc.cmd_array()[0], required: false)
    if _cc_path.found() and host_machine.system() == 'linux' and import('fs').exists('/usr/bin/ld')
      _cc_full = _cc_path.full_path()
      if _conda_prefix == '' or not _cc_full.startswith(_conda_prefix + '/')
        _fuzz_bindir = ['-B/usr/bin/']
      endif
    endif
    foreach _fuzz : [
      ['capi_fuzzer', 'capi_read_fuzzer'],
      ['capi_array_fuzzer', 'capi_array_fuzzer'],
      ['capi_reload_fuzzer', 'capi_reload_fuzzer'],
    ]
      _fuzz_exe = custom_target(
        _fuzz[1],
        input: files('fuzz' / (_fuzz[1] + '.c'), 'fuzz' / 'glibc_isoc23.c'),
        output: _fuzz[1],
        command: _fuzz_cc.cmd_array() + _fuzz_bindir + [
          _fuzz_cc.get_id() == 'clang-cl' ? '/std:c11' : '-std=c11', '-g', '-O1',
          '-I' + meson.project_source_root() / 'include',
          '-I' + meson.current_build_dir() / 'src',
          # The Debian libFuzzer archive is built against the host glibc.
          # A conda gcc links an older libc, which has no __isoc23_strto*.
          # @INPUT@ is the fuzzer source and that shim, so either edit rebuilds.
          '@INPUT@',
          '-o', '@OUTPUT@',
          seams_api_lib.full_path(),
        ] + _fuzz_args + [
          '-Wl,-rpath,' + meson.current_build_dir() / 'src',
        ] + (_omp_rpath != '' ? [_omp_rpath] : []) + [
          # libseams is passed as a file, so the driver does not add the
          # OpenMP runtime it was built with. The harness needs it, and on
          # Darwin that runtime is @rpath/libomp.dylib.
          '-fopenmp',
          '-lstdc++', '-lpthread', '-lm',
        ],
        depends: [seams_api_lib],
        build_by_default: true)
      test(
        _fuzz[0],
        _fuzz_exe,
        env: _fuzz_env,
        args: [
          '-max_total_time=' + get_option('capi_fuzz_seconds').to_string(),
          '-timeout=15', '-rss_limit_mb=4096', '-fork=1', '-ignore_timeouts=1',
        ],
        timeout: get_option('capi_fuzz_seconds') + 180,
        is_parallel: false)
    endforeach
  endif
endif
