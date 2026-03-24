# d-SEAMS Python bindings (nanobind)

py = import('python').find_installation(pure: false)
nanobind_dep = dependency('nanobind')

# --enable-new-dtags is Linux-only (macOS ld does not support it)
_py_link_args = []
if host_machine.system() == 'linux'
  _py_link_args += ['-Wl,--enable-new-dtags']
endif

py.extension_module('_core',
  'bindings.cpp',
  dependencies: [yds_dep, nanobind_dep],
  cpp_args: _args,
  link_args: _py_link_args,
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
