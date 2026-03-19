project('vesin', 'cpp',
  version: '0.3.0',
  default_options: ['cpp_std=c++17'])

vesin_sources = files(
  'vesin/src/vesin.cpp',
  'vesin/src/cpu_cell_list.cpp',
  'vesin/src/shared_libraries.cpp',
  'vesin_cuda_stub.cpp',
)

vesin_inc = include_directories('vesin/include', 'vesin/src')

# CPU cell-list only. GPU path requires vesin's CMake build with gpu-lite
# to generate CUDA kernel code. The stub satisfies linker symbols and
# throws at runtime if CUDA device is requested.
vesin_lib = static_library('vesin',
  vesin_sources,
  include_directories: vesin_inc,
  install: false)

vesin_dep = declare_dependency(
  link_with: vesin_lib,
  include_directories: include_directories('vesin/include'))
