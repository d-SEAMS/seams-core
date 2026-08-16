project('gpulite', 'cpp',
  version: '0.1.0',
  default_options: ['cpp_std=c++17'])

_cc = meson.get_compiler('cpp')
_dl = _cc.find_library('dl', required: host_machine.system() != 'windows')
gpulite_dep = declare_dependency(
  include_directories: include_directories('.'),
  dependencies: [dependency('threads'), _dl])
