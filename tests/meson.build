# d-SEAMS tests (Catch2 v3)

catch2_dep = dependency('catch2-with-main', required: true)

test_sources = {
  'franzblau': 'franzblau-test.cpp',
  'topo_bulk': 'topo_bulk-test.cpp',
  'topo_one_dim': 'topo_one_dim-test.cpp',
  'absor': 'absor-test.cpp',
}

foreach name, src : test_sources
  test(name,
    executable('test_' + name, src,
      dependencies: [yds_dep, catch2_dep],
      cpp_args: _args,
      include_directories: _incdirs),
    workdir: meson.project_source_root() / 'input')
endforeach
