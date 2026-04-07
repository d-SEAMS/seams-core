# d-SEAMS core library (yodaLib)

ydslib_sources = files(
  'absOrientation.cpp',
  'backward.cpp',
  'bond.cpp',
  'bop.cpp',
  'bulkTUM.cpp',
  'cluster.cpp',
  'franzblau.cpp',
  'generic.cpp',
  'mol_sys.cpp',
  'neighbours.cpp',
  'order_parameter.cpp',
  'pntCorrespondence.cpp',
  'rdf2d.cpp',
  'ring.cpp',
  'seams_input.cpp',
  'seams_output.cpp',
  'selection.cpp',
  'shapeMatch.cpp',
  'topo_bulk.cpp',
  'topo_one_dim.cpp',
  'topo_two_dim.cpp',
)

ydslib = library('yodaLib',
  ydslib_sources,
  dependencies: _deps,
  cpp_args: _args,
  include_directories: _incdirs,
  install: true)

_linkto += ydslib

yds_dep = declare_dependency(
  link_with: [ydslib],
  include_directories: _incdirs,
  compile_args: _args,
  dependencies: _deps)
