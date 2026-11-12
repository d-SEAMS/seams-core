project('dseams', 'cpp', version: '2.14.0', default_options: ['cpp_std=c++17'])

py = import('python').find_installation(pure: false)
nanobind_dep = dependency('nanobind')
seams_dep = dependency('seams')

py.extension_module(
  '_core',
  'dseams/_core.cpp',
  dependencies: [nanobind_dep, seams_dep],
  install: true,
  subdir: 'dseams')

py.install_sources(
  ['dseams/__init__.py', 'dseams/__main__.py', 'frames.py'],
  subdir: 'dseams',
  pure: true)
