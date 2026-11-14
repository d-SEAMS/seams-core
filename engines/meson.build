# Engine drivers. Each target is created when its library is found.
# with_lammps and with_gromacs default to auto. A found library whose
# plugin cannot be compiled is an error, not a skipped test.

fs = import('fs')
seams_inc = include_directories('../include')
# seams.h includes the generated seams_version.h.
seams_api_args = ['-I' + meson.project_build_root() / 'src']

lammps_opt = get_option('with_lammps')
if lammps_opt.allowed()
  lammps_dep = dependency('lammps', required: false)
  if not lammps_dep.found()
    lammps_dep = cppc.find_library('lammps', has_headers: ['lammps/library.h'],
                                  required: false)
  endif
  if not lammps_dep.found() and lammps_opt.enabled()
    error('with_lammps is enabled but liblammps was not found')
  endif
  if lammps_dep.found()
    # run_hook.cpp includes "library.h". conda-forge and EESSI keep the
    # headers in include/lammps, one level below the include path.
    lammps_incs = [seams_inc, include_directories('lammps')]
    lammps_header_dir = ''
    lammps_py = find_program('python3', required: false)
    if lammps_py.found()
      lammps_where = run_command(lammps_py, '-c',
        'import os, pathlib, sys\n' +
        'dirs = []\n' +
        'for var in ("CPLUS_INCLUDE_PATH", "CPATH"):\n' +
        '    dirs.extend(pathlib.Path(p) for p in os.environ.get(var, "").split(os.pathsep) if p)\n' +
        'prefix = os.environ.get("CONDA_PREFIX", "")\n' +
        'if prefix:\n' +
        '    dirs.append(pathlib.Path(prefix) / "include")\n' +
        'dirs.extend((pathlib.Path("/usr/include"), pathlib.Path("/usr/local/include")))\n' +
        'for root in dirs:\n' +
        '    header = root / "lammps" / "library.h"\n' +
        '    if header.is_file():\n' +
        '        print(header.parent)\n' +
        '        raise SystemExit(0)\n' +
        'raise SystemExit(1)\n',
        check: false)
      if lammps_where.returncode() == 0
        lammps_header_dir = lammps_where.stdout().strip()
      endif
    endif
    # cpp_args, not include_directories: a pixi env lives inside the
    # source tree, and Meson rejects that absolute include dir.
    lammps_args = lammps_header_dir != '' ? ['-I' + lammps_header_dir] : []
    if not cppc.has_header('library.h', dependencies: lammps_dep,
                           include_directories: lammps_incs,
                           args: lammps_args, required: false)
      error('liblammps was found but library.h is not on the include path')
    endif
    lammps_rpath = ''
    if lammps_dep.type_name() == 'pkgconfig'
      lammps_rpath = lammps_dep.get_variable('libdir', default_value: '')
    endif
    if lammps_rpath == ''
      lammps_lib_py = find_program('python3', required: false)
      if lammps_lib_py.found()
        lammps_lib_where = run_command(lammps_lib_py, '-c',
          'import os, pathlib\n' +
          'prefix = os.environ.get("CONDA_PREFIX", "")\n' +
          'if not prefix:\n' +
          '    raise SystemExit(1)\n' +
          'root = pathlib.Path(prefix)\n' +
          'for folder in (root / "lib", root / "lib64"):\n' +
          '    if (folder / "liblammps.so").is_file() or (folder / "liblammps.so.0").is_file():\n' +
          '        print(folder)\n' +
          '        raise SystemExit(0)\n' +
          'raise SystemExit(1)\n',
          check: false)
        if lammps_lib_where.returncode() == 0 and lammps_lib_where.stdout().strip() != ''
          lammps_rpath = lammps_lib_where.stdout().strip()
        endif
      endif
    endif
    if lammps_rpath == '' and lammps_header_dir != ''
      lammps_prefix = fs.parent(fs.parent(lammps_header_dir))
      foreach sub : ['lib64', 'lib']
        if lammps_rpath == '' and fs.exists(lammps_prefix / sub / 'liblammps.so')
          lammps_rpath = lammps_prefix / sub
        endif
      endforeach
    endif
    run_hook = executable(
      'run_hook',
      files('lammps/run_hook.cpp'),
      dependencies: [lammps_dep],
      link_with: seams_api_lib,
      include_directories: lammps_incs,
      cpp_args: _args + lammps_args + seams_api_args,
      build_rpath: lammps_rpath,
      install_rpath: lammps_rpath,
      install: true)
    test(
      'lammps_hook',
      find_program('skip_missing_lib.sh'),
      args: [run_hook.full_path(), 'in.dseams'],
      workdir: meson.project_source_root() / 'engines/lammps',
      timeout: 120)
    # compute dseams as a `plugin load` module. The LAMMPS class headers
    # include mpi.h, from the MPI library or from the serial STUBS.
    # The host lmp already holds libmpi, so the module takes the headers only.
    lammps_mpi = dependency('mpi', language: 'cpp', required: false)
    lammps_plugin_deps = [lammps_dep]
    if lammps_mpi.found()
      lammps_plugin_deps += lammps_mpi.partial_dependency(compile_args: true, includes: true)
    endif
    lammps_plugin_args = lammps_args
    if get_option('lammps_src_include') != ''
      lammps_plugin_args += ['-I' + get_option('lammps_src_include')]
    endif
    # check_header compiles the include, so a src tree without fmt/ fails here.
    if cppc.check_header('compute.h', dependencies: lammps_plugin_deps,
                         include_directories: lammps_incs, args: lammps_plugin_args,
                         required: false)
      lammps_plugin = shared_module(
        'dseamsplugin',
        files('lammps/compute_dseams.cpp', 'lammps/plugin_dseams.cpp'),
        dependencies: lammps_plugin_deps,
        link_with: seams_api_lib,
        include_directories: lammps_incs,
        cpp_args: _args + lammps_plugin_args + seams_api_args,
        name_prefix: '',
        build_rpath: lammps_rpath,
        install_rpath: lammps_rpath,
        install: true)
    else
      message('LAMMPS class headers are not installed; set lammps_src_include to build compute dseams')
    endif
  endif
endif

gromacs_opt = get_option('with_gromacs')
if gromacs_opt.allowed()
  gromacs_rpath = ''
  gromacs_dep = dependency('libgromacs', required: false)
  if gromacs_dep.found()
    gromacs_rpath = gromacs_dep.get_variable('libdir', default_value: '')
  endif
  if not gromacs_dep.found()
    gromacs_dep = cppc.find_library('gromacs', has_headers: ['gromacs/trajectoryanalysis/analysismodule.h'],
                                   required: false)
  endif
  # conda-forge ships libgromacs under lib.AVX2_256, not lib/.
  if not gromacs_dep.found()
    gromacs_py = find_program('python3', required: false)
    if gromacs_py.found()
      gromacs_where = run_command(gromacs_py, '-c',
        'import os, pathlib, sys\n' +
        'prefix = os.environ.get("CONDA_PREFIX", "")\n' +
        'if not prefix:\n' +
        '    raise SystemExit(1)\n' +
        'root = pathlib.Path(prefix)\n' +
        'names = ("lib.AVX2_256", "lib.AVX_256", "lib")\n' +
        'for name in names:\n' +
        '    lib = root / name / "libgromacs.so"\n' +
        '    if lib.is_file():\n' +
        '        print(lib.parent)\n' +
        '        raise SystemExit(0)\n' +
        'raise SystemExit(1)\n',
        check: false)
      if gromacs_where.returncode() == 0 and gromacs_where.stdout().strip() != ''
        gromacs_rpath = gromacs_where.stdout().strip()
        gromacs_dep = cppc.find_library('gromacs',
          dirs: [gromacs_rpath],
          has_headers: ['gromacs/trajectoryanalysis/analysismodule.h'],
          required: false)
      endif
    endif
  endif
  if not gromacs_dep.found() and gromacs_opt.enabled()
    error('with_gromacs is enabled but libgromacs was not found')
  endif
  if gromacs_dep.found() and gromacs_rpath == ''
    gromacs_py = find_program('python3', required: false)
    if gromacs_py.found()
      gromacs_where = run_command(gromacs_py, '-c',
        'import os, pathlib\n' +
        'prefix = os.environ.get("CONDA_PREFIX", "")\n' +
        'if not prefix:\n' +
        '    raise SystemExit(1)\n' +
        'root = pathlib.Path(prefix)\n' +
        'names = ("lib.AVX2_256", "lib.AVX_256", "lib")\n' +
        'for name in names:\n' +
        '    lib = root / name / "libgromacs.so"\n' +
        '    if lib.is_file():\n' +
        '        print(lib.parent)\n' +
        '        raise SystemExit(0)\n' +
        'raise SystemExit(1)\n',
        check: false)
      if gromacs_where.returncode() == 0 and gromacs_where.stdout().strip() != ''
        gromacs_rpath = gromacs_where.stdout().strip()
      endif
    endif
  endif
  gromacs_args = []
  if gromacs_dep.found() and get_option('gromacs_legacy_include') != ''
    gromacs_args += ['-I' + get_option('gromacs_legacy_include'),
      '-I' + (meson.current_source_dir() / 'gromacs' / 'compat')]
  endif
  gromacs_header = gromacs_dep.found() and cppc.has_header(
    'gromacs/trajectoryanalysis/analysismodule.h', dependencies: gromacs_dep,
    args: gromacs_args, required: false)
  if gromacs_dep.found() and not gromacs_header
    if gromacs_opt.enabled()
      error('libgromacs was found without its trajectory-analysis headers; set gromacs_legacy_include')
    endif
    message('libgromacs has no trajectory-analysis headers; set gromacs_legacy_include to build gmx-seams')
  endif
  if gromacs_header
    gmx_seams = executable(
      'gmx-seams',
      files('gromacs/seams_ice.cpp'),
      dependencies: [gromacs_dep],
      link_with: seams_api_lib,
      include_directories: [seams_inc, include_directories('.')],
      cpp_args: _args + seams_api_args + gromacs_args,
      # GROMACS 2024 headers call ssize() unqualified, ambiguous with std::ssize.
      override_options: ['cpp_std=c++17'],
      build_rpath: gromacs_rpath,
      install_rpath: gromacs_rpath,
      install: true)
    test(
      'gromacs_mdrun',
      find_program('gromacs/run_mdrun.sh'),
      args: [gmx_seams.full_path(), fs.parent(seams_api_lib.full_path())],
      timeout: 180)
  endif
endif

openmm_opt = get_option('with_openmm')
if openmm_opt.allowed()
  openmm_dep = dependency('OpenMM', required: false)
  if not openmm_dep.found()
    openmm_dep = cppc.find_library('OpenMM', has_headers: ['OpenMM.h'], required: false)
  endif
  if not openmm_dep.found() and openmm_opt.enabled()
    error('with_openmm is enabled but libOpenMM was not found')
  endif
  if openmm_dep.found()
    # The test runs the driver without LD_LIBRARY_PATH.
    openmm_rpath = ''
    if openmm_dep.type_name() == 'pkgconfig'
      openmm_rpath = openmm_dep.get_variable('libdir', default_value: '')
    endif
    if openmm_rpath == ''
      openmm_lib_py = find_program('python3', required: false)
      if openmm_lib_py.found()
        openmm_lib_where = run_command(openmm_lib_py, '-c',
          'import os, pathlib\n' +
          'dirs = []\n' +
          'for var in ("CONDA_PREFIX", "OPENMM_DIR"):\n' +
          '    if os.environ.get(var):\n' +
          '        dirs += [pathlib.Path(os.environ[var]) / "lib", pathlib.Path(os.environ[var]) / "lib64"]\n' +
          'for var in ("LIBRARY_PATH", "LD_RUN_PATH"):\n' +
          '    dirs.extend(pathlib.Path(p) for p in os.environ.get(var, "").split(os.pathsep) if p)\n' +
          'for root in dirs:\n' +
          '    if (root / "libOpenMM.so").is_file():\n' +
          '        print(root)\n' +
          '        raise SystemExit(0)\n' +
          'raise SystemExit(1)\n',
          check: false)
        if openmm_lib_where.returncode() == 0 and openmm_lib_where.stdout().strip() != ''
          openmm_rpath = openmm_lib_where.stdout().strip()
        endif
      endif
    endif
    openmm_plugin = executable(
      'seams_openmm',
      files('openmm/plugin.cpp'),
      dependencies: [openmm_dep],
      link_with: seams_api_lib,
      include_directories: seams_inc,
      cpp_args: _args + seams_api_args,
      build_rpath: openmm_rpath,
      install_rpath: openmm_rpath,
      install: true)
    test(
      'openmm_plugin',
      find_program('skip_missing_lib.sh'),
      args: [openmm_plugin.full_path()],
      timeout: 120)
  endif
endif

plumed_opt = get_option('with_plumed')
if plumed_opt.allowed()
  plumed_dep = dependency('plumed', required: false)
  if not plumed_dep.found()
    plumed_dep = cppc.find_library('plumed', has_headers: ['plumed/wrapper/Plumed.h'],
                                  required: false)
  endif
  # The collective variable is libSeamsPlumed.so. Installed headers are
  # plumed/core/Colvar.h under the prefix include directory.
  plumed_prog = find_program('plumed', required: false)
  if plumed_dep.found() or plumed_prog.found()
    plumed_incs = [seams_inc, include_directories('.')]
    plumed_probe_deps = plumed_dep.found() ? [plumed_dep] : []
    # libplumed.so is the loader. The collective-variable symbols live in
    # libplumedKernel.so.
    # The kernel must be the one that matches this PLUMED. A conda prefix
    # earlier on the linker path can supply an older libplumedKernel.
    plumed_libdir = ''
    if plumed_dep.found() and plumed_dep.type_name() == 'pkgconfig'
      plumed_libdir = plumed_dep.get_variable('libdir', default_value: '')
    endif
    if plumed_libdir == '' and plumed_prog.found()
      plumed_libdir = fs.parent(fs.parent(plumed_prog.full_path())) / 'lib'
    endif
    plumed_kernel = cppc.find_library('plumedKernel',
      dirs: plumed_libdir != '' ? [plumed_libdir] : [],
      required: false)
    if plumed_kernel.found()
      plumed_probe_deps += plumed_kernel
    endif
    plumed_header = cppc.has_header('plumed/core/Colvar.h',
                                    dependencies: plumed_probe_deps,
                                    required: false)
    if not plumed_header and plumed_prog.found()
      info = run_command(plumed_prog, 'info', '--include-dir', check: false)
      if info.returncode() == 0 and info.stdout().strip() != ''
        plumed_incs += include_directories(info.stdout().strip())
        plumed_header = cppc.has_header('plumed/core/Colvar.h',
                                        include_directories: plumed_incs,
                                        dependencies: plumed_probe_deps,
                                        required: false)
      endif
    endif
    if not plumed_header
      error('PLUMED was found but plumed/core/Colvar.h is not usable, so libSeamsPlumed.so cannot be built')
    endif
    # CPATH is searched before an -isystem copy of the same directory, so a
    # second PLUMED prefix can supply the headers while this library is linked.
    # A distinct -I of the prefix that owns the library is searched first.
    plumed_incdir = ''
    if plumed_dep.found() and plumed_dep.type_name() == 'pkgconfig'
      plumed_incdir = plumed_dep.get_variable('includedir', default_value: '')
    endif
    if plumed_incdir == '' and plumed_prog.found()
      info = run_command(plumed_prog, 'info', '--include-dir', check: false)
      if info.returncode() == 0
        plumed_incdir = info.stdout().strip()
      endif
    endif
    plumed_cpp = []
    if plumed_incdir != ''
      # A symlink of the prefix include is the same directory as -isystem,
      # and GCC drops that -I. A directory that only holds the plumed tree
      # stays distinct and is searched before CPATH.
      plumed_link = meson.project_build_root() / 'plumed-include'
      run_command('mkdir', '-p', plumed_link, check: false)
      run_command('ln', '-sfn', plumed_incdir / 'plumed', plumed_link / 'plumed', check: false)
      plumed_cpp = ['-I' + plumed_link]
    endif
    seams_plumed = shared_library(
      'SeamsPlumed',
      files('plumed/SeamsColvar.cpp', 'plumed/plumed_components.cpp'),
      dependencies: plumed_probe_deps,
      link_with: seams_api_lib,
      include_directories: plumed_incs,
      cpp_args: _args + seams_api_args + plumed_cpp,
      name_prefix: 'lib',
      # PLUMED's LOAD dlopens this; without an rpath libplumedKernel.so
      # resolves only through LD_LIBRARY_PATH.
      build_rpath: plumed_libdir,
      install_rpath: plumed_libdir,
      install: true)
  elif plumed_opt.enabled()
    error('with_plumed is enabled but PLUMED was not found')
  endif
endif

triclinic_box_test = executable(
  'test_triclinic_box',
  files('../tests/test_triclinic_box.cpp'),
  cpp_args: _args)
test('triclinic_box', triclinic_box_test, timeout: 30)

python3_engines = find_program('python3', required: false)
if python3_engines.found() and get_option('with_tests')
  openmm_env = environment()
  openmm_env.set('DSEAMS_BUILD', meson.project_build_root())
  openmm_env.set('DSEAMS_SRC', meson.project_source_root())
  openmm_env.set('PYTHONPATH', meson.project_build_root() / 'python-stage')
  test(
    'openmm_reporter',
    python3_engines,
    args: [meson.project_source_root() / 'engines/openmm/run_reporter.py'],
    env: openmm_env,
    timeout: 180)
  test(
    'plumed_driver',
    python3_engines,
    args: [
      meson.project_source_root() / 'engines/plumed/run_driver.py',
      meson.project_build_root(),
    ],
    timeout: 180)
  test(
    'insitu_observe',
    python3_engines,
    args: [
      meson.project_source_root() / 'engines/insitu_observe.py',
      meson.project_build_root(),
    ],
    timeout: 120)
endif
