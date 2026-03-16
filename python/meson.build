# d-SEAMS Python bindings (nanobind)

py = import('python').find_installation(pure: false)
nanobind_dep = dependency('nanobind')

py.extension_module('_core',
  'bindings.cpp',
  dependencies: [yds_dep, nanobind_dep],
  cpp_args: _args,
  install: true,
  subdir: 'pydseamslib')

py.install_sources(
  ['pydseamslib/__init__.py'],
  pure: false,
  subdir: 'pydseamslib')

py.install_sources(
  ['pydseamslib/adapters/__init__.py',
   'pydseamslib/adapters/_ase.py'],
  pure: false,
  subdir: 'pydseamslib/adapters')
