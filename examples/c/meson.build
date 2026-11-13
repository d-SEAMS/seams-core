project('seams-c-quickstart', 'c', default_options: ['c_std=c11'])

seams_dep = dependency('seams')
m_dep = meson.get_compiler('c').find_library('m', required: false)
executable('quickstart', 'quickstart.c', dependencies: [seams_dep, m_dep])
