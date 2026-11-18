project('vesin', 'cpp',
  version: '0.5.7',
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
# The parent coverage build instruments this static library. Neighbour
# lists run on several threads, and a plain arc update then goes negative.
vesin_cpp_args = []
if get_option('b_coverage')
  vesin_cpp_args += ['-fprofile-update=atomic']
endif
vesin_lib = static_library('vesin',
  vesin_sources,
  include_directories: vesin_inc,
  cpp_args: vesin_cpp_args,
  # Upstream sources; a parent -Dwerror=true must not fail on their warnings.
  override_options: ['werror=false'],
  install: false)

vesin_dep = declare_dependency(
  link_with: vesin_lib,
  include_directories: include_directories('vesin/include'))
