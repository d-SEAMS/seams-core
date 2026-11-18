# Generated beside this directory in the build tree. A file under the
# build tree's tests/ becomes a second "tests" directory once genhtml
# strips the build prefix and the source prefix, and the report stops.
_cuda_kernels = custom_target('gpu_batch_kernels',
  input: meson.project_source_root() / 'src' / 'gpu_batch.cpp',
  output: 'gpu_batch_kernels.inc',
  command: [find_program('python3'),
            meson.project_source_root() / 'tests' / 'extract_cuda_kernels.py',
            '@INPUT@', '@OUTPUT@'])
_gpu_host_kernels_inc = include_directories('.')
