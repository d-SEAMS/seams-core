# structureFactor

[![built with nix](https://builtwithnix.org/badge.svg)](https://builtwithnix.org)

Check our docs build status [here](https://travis-ci.org/d-SEAMS/seams-core).
The docs themselves are [here](https://d-seams.github.io/seams-core/) and
development is ongoing [on GitHub](https://github.com/d-SEAMS/seams-core).

This C++14 program reads in lammps trajectory files
or XYZ files to calculate RDF, in-plane RDF and the structure factor.

\brief Post-processing code for RDF, in-plane RDF, structure factor

# Compilation

## Dependency Management

Unless you are sure that all the libraries necessary are available on your
system, it is expedient to ensure `conan` is available for managing
dependencies. Since `conan` is essentially a python package, it is handled
by[poetry](https://github.com/sdispater/poetry).

### Python

More recent versions use [poetry](https://github.com/sdispater/poetry) which as
a non-system python dependant installation.

```{bash}
# Get poetry (once)
curl -sSL https://raw.githubusercontent.com/sdispater/poetry/master/get-poetry.py | python
# Install packages
poetry install
```

### Lua

TODO: Move to conan

Lua v5.3 is used for the scripting engine. It needs to be installed via the
operating system's normal packaging system for now. If possible, install a
version compiled with `c++`, not `c`.

```{bash}
# Ubuntu and derivatives
sudo apt install lua5.3 liblua5.3
# ArchLinux
sudo pacman -S lua
```

### Lua Modules

Since a major portion of the frontend is in `lua`, the following modules are
required.[LuaRocks](https://luarocks.org/) is the recommended package manager
and they are to be installed as root.

```sh
# For cross-OS filesystem operations
sudo luarocks install luafilesystem
```

## Build Things

To generate the executable _yodaStruct_ in `bin/` use:

```{bash}
# Use the bundled conan
poetry shell
# Always prefer an out of tree build
mkdir build
cd build

#
# OPTIONAL STEP (TODO: Use CMakeLists.txt to figure this out)
#
# We prefer clang, so if you have it, use it
export CXX=/usr/bin/clang++
export CC=/usr/bin/clang
#
# END OPTIONAL
#

# If you forget the flag you will build the Debug version
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

# Nix Usage

We use a deterministic build system to generate both bug reports and uniform
usage statistics.

## Build

Since this project is built with `nix`, we can simply do the following from the
root directory:

```sh
# This will take a long time the first time since it builds sharkML
nix-build .
# Install into your path
nix-env -if .
# Use anywhere
cd lua_inputs/
yodaStruct -c config.yml
# Use with lua modules
nix-shell --run 'bash' --pure
```

### Caveats

Though the build itself is guaranteed to be reproducible as the `nixpkgs` are
also pinned to a particular commit, the `luarocks` dependencies are still local,
since they are determined at runtime. This means, for example, to use the sample
file, you need to ensure you have the `luarocks` modules installled in your
system.

The above caveats are not relevant when you run it in the shell environment
defined by `shell.nix`

#### Reproducible Lua

For reproducing `lua` we use [luas](https://github.com/limadm/luas). Note that
this is still an imperfect method and the best way to run this is via the
`nix-shell --run 'bash' --pure` environment.

```sh
luas init 5.2.4
luas use 5.2.4
luarocks install luafilesystem
```

## Development

We can simply use the `nix` environment:

```sh
# From the project root
nix-shell
# Sanitize and fix the shell
stty sane
export TERM="xterm-256color"
```

# Running

To run the sample inputs, simply move the binary to the project root, or to a
directory where `input/` is a child directory.

```{bash}
# Assuming you are in the build directory
# Check help with -h
# --script and --file are optional now
./yodaStruct --script ../lua_inputs/transition_diff.lua -f ../lua_inputs/parameter.txt -c ../lua_inputs/config.yml
```

# Details

- For example, if you want to find the in-plane RDF in the
  XY plane, you should use the Rdf2D class. The equation used for 2D-RDF for the \f$n^{th}\f$ layer is:
  \f[
  g^n(r) = \frac{1}{(\rho^n)^2 A \delta z} \Sigma*{i \neq j} \delta(r - r*{ij}) \left[ \Theta\left( \frac{\delta z}{2}-|z_i-z^2| \right) \times \Theta\left( \frac{\delta z}{2}-|z_j-z^n| \right) \right]
  \f]
  For detailed instructions, see the Rdf2D class documentation

# Developer Documentation

TODO: Move this to some other location.

For updates to any of the **bundled** `external libraries` change the commit number and use:

```{bash}
$ cd src/external
# Sol2
 wget https://raw.githubusercontent.com/ThePhD/sol2/develop/single/sol/sol_forward.hpp
 wget https://raw.githubusercontent.com/ThePhD/sol2/develop/single/sol/sol.hpp
# cxxopts
 wget https://raw.githubusercontent.com/jarro2783/cxxopts/master/include/cxxopts.hpp
```

## Leaks and performance

While testing for leaks, use `clang` (for
[AddressSanitizer](https://github.com/google/sanitizers/wiki/AddressSanitizer)
and [LeakSanitizer](https://github.com/google/sanitizers/wiki/AddressSanitizerLeakSanitizer))
and the following:

```{bash}
export CXX=/usr/bin/clang++ && export CC=/usr/bin/clang
cmake .. -DCMAKE_CXX_FLAGS="-pg -fsanitize=address " -DCMAKE_EXE_LINKER_FLAGS=-pg -DCMAKE_SHARED_LINKER_FLAGS=-pg
```

# Contributing

Please ensure that all contributions are formatted according to the
[clang-format](./clang-format) configuration file.

Specifically, consider using the following:

-[Sublime Plugin](https://github.com/rosshemsley/SublimeClangFormat) for users
of Sublime Text

- [format-all](https://github.com/lassik/emacs-format-all-the-code) for Emacs
- [vim-clang-format](https://github.com/rhysd/vim-clang-format) for Vim
- Visual Studio: http://llvm.org/builds/, or use the [integrated support in Visual Studio 2017](https://blogs.msdn.microsoft.com/vcblog/2018/03/13/clangformat-support-in-visual-studio-2017-15-7-preview-1/)
- Xcode: https://github.com/travisjeffery/ClangFormat-Xcode

Where some of the above suggestions are derived from [this depreciated githook](https://github.com/andrewseidl/githook-clang-format).

# Acknowledgements

The following tools are used in this project:

- [CMake](https://cmake.org/) for compilation ([cmake-init](https://github.com/cginternals/cmake-init) was used as a reference)
- [Clang](https://clang.llvm.org/) because it is more descriptive with better tools
- [Conan](https://conan.io/) and [https://pipenv.readthedocs.io/en/latest/](pipenv) for dependency management
- [Doxygen](https://www.doxygen.org) for the developer API
- [clang-format](https://clang.llvm.org/docs/ClangFormat.html) for code formatting
- [lua](https://www.lua.org) for the scripting engine
- [yaml](http://yaml.org/) for the configuration

## Third Party Libraries

The libraries used are:

- [backward-cpp](https://github.com/bombela/backward-cpp) for better stacktraces without `gdb`
- [cxxopts](https://github.com/jarro2783/cxxopts) for parsing command line options
- [rang](https://github.com/agauniyal/rang) for terminal styles (ANSI)
- [sol2](https://github.com/ThePhD/sol2) for interfacing with lua
- [yaml-cpp](https://github.com/jbeder/yaml-cpp) for working with `yaml`
- [fmt](https://github.com/fmtlib/fmt) for safe and fast formatting
- [Linear Algebra PACKage (LAPACK)](http://www.netlib.org/lapack/)
- [Basic Linear Algebra Subprograms (BLAS)](http://www.netlib.org/blas/)
- [Boost Geometry](https://www.boost.org/doc/libs/1_68_0/libs/geometry/doc/html/index.html) for working with different coordinates
- [Boost Math](https://www.boost.org/doc/libs/?view=category_math) for spherical harmonics
- [Blaze](https://bitbucket.org/blaze-lib/blaze/) for very fast modern linear algebra
- [nanoflann](https://github.com/jlblancoc/nanoflann) to calculate nearest neighbors
