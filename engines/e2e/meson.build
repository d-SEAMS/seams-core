# Engine end-to-end checks. Each runs a short ice/liquid job in the real
# engine, labels it in situ, and compares every label of every dumped frame
# with `seams chill-plus` on the same frames. A missing engine is a SKIP
# (exit 77). Engine executables come from PATH or LMP, GMX, PLUMED_KERNEL.

e2e_py = find_program('python3', required: false)
if e2e_py.found()
  e2e_src = meson.project_source_root() / 'engines'
  e2e_build = meson.project_build_root()
  e2e_env = environment()
  e2e_env.set('SEAMS_CLI', seams_exe.full_path())
  e2e_deps = [seams_exe, seams_api_lib]
  foreach name : ['lammps_plugin', 'gmx_seams', 'seams_plumed']
    if is_variable(name)
      e2e_deps += get_variable(name)
    endif
  endforeach
  test('lammps_e2e', e2e_py,
    args: [e2e_src / 'lammps' / 'e2e_mw.py', e2e_build],
    env: e2e_env, depends: e2e_deps, suite: 'engines-e2e', timeout: 900)
  # Regression: on two ranks each rank used to label only its own atoms.
  test('lammps_e2e_np2', e2e_py,
    args: [e2e_src / 'lammps' / 'e2e_mw.py', e2e_build, '--ranks', '2'],
    env: e2e_env, depends: e2e_deps, suite: 'engines-e2e', timeout: 900)
  # Regression: compute dseams had no type or mol keyword, so a frame with
  # hydrogens or ions put them in the CHILL+ graph.
  test('lammps_select', e2e_py,
    args: [e2e_src / 'lammps' / 'e2e_select.py', e2e_build],
    env: e2e_env, depends: e2e_deps, suite: 'engines-e2e', timeout: 300)
  test('lammps_select_np2', e2e_py,
    args: [e2e_src / 'lammps' / 'e2e_select.py', e2e_build, '--ranks', '2'],
    env: e2e_env, depends: e2e_deps, suite: 'engines-e2e', timeout: 300)
  test('plumed_e2e', e2e_py,
    args: [e2e_src / 'plumed' / 'e2e_mw.py', e2e_build],
    env: e2e_env, depends: e2e_deps, suite: 'engines-e2e', timeout: 900)
  # Regression: SEAMS assumed PLUMED lengths were nm.
  test('plumed_units', e2e_py,
    args: [e2e_src / 'plumed' / 'run_units.py', e2e_build],
    depends: e2e_deps, suite: 'engines-e2e', timeout: 120)
  test('gromacs_e2e', e2e_py,
    args: [e2e_src / 'gromacs' / 'e2e_tip4pice.py', e2e_build],
    env: e2e_env, depends: e2e_deps, suite: 'engines-e2e', timeout: 900)
  test('gromacs_plumed_e2e', e2e_py,
    args: [e2e_src / 'gromacs' / 'e2e_tip4pice.py', e2e_build, '--plumed'],
    env: e2e_env, depends: e2e_deps, suite: 'engines-e2e', timeout: 900)
  test('openmm_e2e', e2e_py,
    args: [e2e_src / 'openmm' / 'e2e_tip4pice.py', e2e_build],
    env: e2e_env, depends: e2e_deps, suite: 'engines-e2e', timeout: 900)
endif
