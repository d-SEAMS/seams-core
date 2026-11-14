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
    # run_hook.cpp includes "library.h". conda-forge keeps that header in
    # include/lammps, next to the pkg-config includedir.
    lammps_incs = [seams_inc, include_directories('lammps')]
    lammps_header_arg = ''
    if not cppc.has_header('library.h', dependencies: lammps_dep, required: false)
      lammps_py = find_program('python3', required: false)
      if lammps_py.found()
        lammps_where = run_command(lammps_py, '-c',
          'import os, pathlib, sys\n' +
          'roots = []\n' +
          'prefix = os.environ.get("CONDA_PREFIX", "")\n' +
          'if prefix:\n' +
          '    roots.append(pathlib.Path(prefix))\n' +
          'roots.extend((pathlib.Path("/usr"), pathlib.Path("/usr/local")))\n' +
          'for root in roots:\n' +
          '    header = root / "include" / "lammps" / "library.h"\n' +
          '    if header.is_file():\n' +
          '        print(header.parent)\n' +
          '        raise SystemExit(0)\n' +
          'raise SystemExit(1)\n',
          check: false)
        if lammps_where.returncode() == 0 and lammps_where.stdout().strip() != ''
          # cpp_args, not include_directories: a pixi env lives inside the
          # source tree, and Meson rejects that absolute include dir.
          lammps_header_arg = '-I' + lammps_where.stdout().strip()
        endif
      endif
    endif
    lammps_args = lammps_header_arg != '' ? [lammps_header_arg] : []
    if not cppc.has_header('library.h', dependencies: lammps_dep,
                           include_directories: lammps_incs,
                           args: lammps_args, required: false)
      error('liblammps was found but library.h is not on the include path')
    endif
    run_hook = executable(
      'run_hook',
      files('lammps/run_hook.cpp'),
      dependencies: [lammps_dep],
      link_with: seams_api_lib,
      include_directories: lammps_incs,
      cpp_args: _args + lammps_args + seams_api_args,
      install: true)
    test(
      'lammps_hook',
      run_hook,
      args: ['in.dseams'],
      workdir: meson.project_source_root() / 'engines/lammps',
      timeout: 120)
  endif
endif

gromacs_opt = get_option('with_gromacs')
if gromacs_opt.allowed()
  gromacs_rpath = ''
  # System include, so the -Werror conversions of a clang build stop at
  # this project's code and not in the GROMACS headers.
  gromacs_dep = dependency('libgromacs', required: false, include_type: 'system')
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
  if gromacs_dep.found()
    gmx_seams = executable(
      'gmx-seams',
      files('gromacs/seams_ice.cpp'),
      dependencies: [gromacs_dep],
      link_with: seams_api_lib,
      include_directories: [seams_inc, include_directories('.')],
      cpp_args: _args + seams_api_args,
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
    openmm_plugin = executable(
      'seams_openmm',
      files('openmm/plugin.cpp'),
      dependencies: [openmm_dep],
      link_with: seams_api_lib,
      include_directories: seams_inc,
      cpp_args: _args + seams_api_args,
      install: true)
    test('openmm_plugin', openmm_plugin, timeout: 120)
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
    seams_plumed = shared_library(
      'SeamsPlumed',
      files('plumed/SeamsColvar.cpp'),
      dependencies: plumed_probe_deps,
      link_with: seams_api_lib,
      include_directories: plumed_incs,
      cpp_args: _args + seams_api_args,
      name_prefix: 'lib',
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
