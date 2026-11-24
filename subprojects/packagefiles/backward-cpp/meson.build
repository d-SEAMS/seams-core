project(
  'backward-cpp',
  'cpp',
  version: '1.6',
  meson_version: '>= 0.47',
)

cpp = meson.get_compiler('cpp')

bfd = cpp.find_library(
  'bfd',
  required: get_option('bfd'),
)
dw = cpp.find_library(
  'dw',
  required: get_option('dw'),
)
libunwind = dependency(
  'libunwind',
  required: get_option('libunwind'),
)
dwarf = dependency(
  'libdwarf',
  required: get_option('dwarf'),
)

deps = []
args = []
# stack walking
if libunwind.found()
  deps += libunwind
  args += ['-DBACKWARD_HAS_LIBUNWIND=1']
endif

# stack detail extraction
if dw.found()
  deps += dw
  args += ['-DBACKWARD_HAS_DW=1']
elif dwarf.found()
  deps += dwarf
  args += ['-DBACKWARD_HAS_DWARF=1']
elif bfd.found()
  deps += bfd
  args += ['-DBACKWARD_HAS_BFD=1']
  if host_machine.system() != 'windows'
    if meson.version().version_compare('>= 0.62')
      deps += dependency('dl')
    else
      deps += cpp.find_library('dl')
    endif
  endif
endif

install_headers(
  'backward.hpp',
  subdir: 'backward',
)

if host_machine.system() == 'windows'
  deps += cpp.find_library('dbghelp')
  backward = static_library(
    'backward-cpp',
    'backward.cpp',
    dependencies: deps,
    install: true,
    cpp_args: args,
  )
else
  backward = library(
    'backward-cpp',
    'backward.cpp',
    dependencies: deps,
    install: true,
    cpp_args: args,
  )
endif

backward_interface_dep = declare_dependency(
  include_directories: include_directories('.'),
  dependencies: deps,
  compile_args: args,
)

backward_dep = declare_dependency(
  dependencies: backward_interface_dep,
  link_with: backward,
)

# This file is a local overlay of WrapDB's backward-cpp_1.6-6 packaging
# (subprojects/packagefiles/backward-cpp, applied through this wrap's
# patch_directory key), not the stock WrapDB patch. It differs from
# WrapDB's own meson.build only by dropping the trailing `tests = [...]`
# / `foreach test : tests ... test(...)` blocks that build
# testexe_test/testexe_stacktrace/testexe_rectrace/testexe_select_signals/
# testexe_suicide as default ninja targets: those are backward-cpp's own
# upstream test suite (intentionally-crashing demo code), registered with
# no `if not meson.is_subproject()` guard, so a consuming project's own
# `ninja`/`-Dwerror=true`/`-Dwith_tests=true` build was building and
# werror-failing them too. This project consumes backward_interface_dep
# only (never backward_dep's compiled `library()`, see src/meson.build),
# so dropping those targets changes nothing this project links against.
