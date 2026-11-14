# Nanobind module over libseams, plus the pure-Python package.
# The caller sets dseams_api_dep: the in-tree C API, or an installed
# libseams when system_libseams is on.

# A development build takes the python3 the tests invoke; Meson's own
# interpreter is not that program when a conda env is on PATH. A
# packaging build (no tests) takes the interpreter the build frontend
# configured, which is the one the wheel is for.
if get_option('with_tests')
  py_mod = import('python').find_installation('python3', pure: false, required: false)
else
  py_mod = import('python').find_installation(pure: false, required: false)
endif

# nanobind's stable ABI starts at 3.12. meson-python only tags the wheel
# abi3 when python.allow_limited_api is true, and that is off by default
# so any interpreter can build the sdist.
dseams_limited_api = ''
if py_mod.found() and get_option('python.allow_limited_api')
  if py_mod.language_version().version_compare('>=3.12') and py_mod.get_variable('Py_GIL_DISABLED', 0) == 0
    dseams_limited_api = '3.12'
  endif
endif

nanobind_dep = dependency('nanobind', required: false, allow_fallback: true,
                          default_options: ['limited_api=' + dseams_limited_api])
python_deps = []
if py_mod.found() and nanobind_dep.found()
  _ext_kwargs = {}
  if dseams_limited_api != ''
    _ext_kwargs += {'limited_api': dseams_limited_api}
  endif
  dseams_core = py_mod.extension_module(
    '_core',
    files('dseams/_core.cpp'),
    dependencies: [nanobind_dep, dseams_api_dep],
    install: true,
    subdir: 'dseams',
    kwargs: _ext_kwargs)
  py_mod.install_sources(
    files(
      'dseams/__init__.py',
      'dseams/__main__.py',
      'dseams/userframe.py',
      'frames.py',
    ),
    subdir: 'dseams',
    pure: false)
  # Documented Frame import and the 2.11 package name, for pydseamslib users.
  py_mod.install_sources(
    files('pydseams/__init__.py', 'pydseams/band.py'),
    subdir: 'pydseams',
    pure: false)
  py_mod.install_sources(
    files('pydseamslib/__init__.py'),
    subdir: 'pydseamslib',
    pure: false)
  python_deps += dseams_core
endif
