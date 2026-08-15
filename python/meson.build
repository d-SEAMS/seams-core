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
  # Resolve the freshly built libyodaLib over a stale installed copy: the
  # conda LDFLAGS rpath to the env prefix otherwise wins
  build_rpath: meson.project_build_root() / 'src',
  install: true,
  subdir: 'pydseamslib',
  # CPython 3.12+ stable ABI. One wheel per platform; not valid on
  # free-threaded interpreters (no limited ABI there).
  limited_api: '3.12')

py.install_sources(
  ['pydseamslib/__init__.py',
   'pydseamslib/frame.py',
   'pydseamslib/aseio.py'],
  pure: false,
  subdir: 'pydseamslib')

py.install_sources(
  ['pydseamslib/adapters/__init__.py',
   'pydseamslib/adapters/_ase.py'],
  pure: false,
  subdir: 'pydseamslib/adapters')
