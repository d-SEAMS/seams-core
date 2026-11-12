# Engine drivers. Each target is created only when its library is found.
# with_lammps and with_gromacs default to disabled so a plain build does
# not look up Ubuntu library paths. The meson C++ compiler links them.

fs = import('fs')
seams_inc = include_directories('../include')

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
    run_hook = executable(
      'run_hook',
      files('lammps/run_hook.cpp'),
      dependencies: [lammps_dep],
      link_with: seams_api_lib,
      include_directories: [seams_inc, include_directories('lammps')],
      cpp_args: _args,
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
  gromacs_dep = dependency('libgromacs', required: false)
  if not gromacs_dep.found()
    gromacs_dep = cppc.find_library('gromacs', has_headers: ['gromacs/trajectoryanalysis/analysismodule.h'],
                                   required: false)
  endif
  if not gromacs_dep.found() and gromacs_opt.enabled()
    error('with_gromacs is enabled but libgromacs was not found')
  endif
  if gromacs_dep.found()
    gmx_seams = executable(
      'gmx-seams',
      files('gromacs/seams_ice.cpp'),
      dependencies: [gromacs_dep],
      link_with: seams_api_lib,
      include_directories: [seams_inc, include_directories('.')],
      cpp_args: _args,
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
      cpp_args: _args,
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
  # The collective variable is a PLUMED plugin (libSeamsPlumed.so). The
  # headers live next to the plumed executable when pkg-config is absent.
  plumed_prog = find_program('plumed', required: false)
  if plumed_dep.found() or plumed_prog.found()
    # Compiled when the compiler can see core/Colvar.h. A missing header
    # leaves the shared module unbuilt; the driver test then skips.
    plumed_header = cppc.has_header('core/Colvar.h', required: false)
    if plumed_header
      seams_plumed = shared_library(
        'SeamsPlumed',
        files('plumed/SeamsColvar.cpp'),
        dependencies: plumed_dep.found() ? [plumed_dep] : [],
        link_with: seams_api_lib,
        include_directories: [seams_inc, include_directories('.')],
        cpp_args: _args,
        name_prefix: 'lib',
        install: true)
    endif
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
  test(
    'openmm_reporter',
    python3_engines,
    args: [meson.project_source_root() / 'engines/openmm/run_reporter.py'],
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
