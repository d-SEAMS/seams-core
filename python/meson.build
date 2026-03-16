# d-SEAMS Python bindings (nanobind)

py = import('python').find_installation(pure: false)
nanobind_dep = dependency('nanobind')

py.extension_module('_core',
  'bindings.cpp',
  dependencies: [yds_dep, nanobind_dep],
  cpp_args: _args,
  link_args: ['-Wl,--enable-new-dtags'],
  install: true,
  subdir: 'pydseamslib',
  limited_api: '3.12')

py.install_sources(
  ['pydseamslib/__init__.py'],
  pure: false,
  subdir: 'pydseamslib')

py.install_sources(
  ['pydseamslib/adapters/__init__.py',
   'pydseamslib/adapters/_ase.py'],
  pure: false,
  subdir: 'pydseamslib/adapters')
