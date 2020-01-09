# d-SEAMS

**Deferred Structural Elucidation Analysis for Molecular Simulations**

[![Docs Status](https://travis-ci.org/d-SEAMS/seams-core.svg?branch=master)](https://travis-ci.org/d-SEAMS/seams-core)
[![built with nix](https://builtwithnix.org/badge.svg)](https://builtwithnix.org)

Check our docs build status [here](https://travis-ci.org/d-SEAMS/seams-core).
The docs themselves are [here](https://d-seams.github.io/seams-core/) and
development is ongoing [on GitHub](https://github.com/d-SEAMS/seams-core). We
also have [a Zenodo community](https://zenodo.org/communities/d-seams/) for user-contributions like reviews, testimonials
and tutorials.

\brief The C++ core of d-SEAMS, a molecular dynamics trajectory analysis engine.

\note The <a href="pages.html">wiki</a> describes the examples and how to obtain
the data-sets (trajectories) <a
href="https://figshare.com/projects/d-SEAMS_Datasets/73545">from figshare</a>.

\alert **If** you are unwilling to use the `nix` build system, then **please note** that you must manage the dependencies MANUALLY, including the compiler versions

# Citation

This software is being actively developed and written up. If you do use it in an
academic capacity, for now please cite [the following preprint](https://arxiv.org/abs/1909.09830):

    Goswami, R.; Goswami, A.; Singh, J. K. (2019). "d-SEAMS: Deferred Structural Elucidation Analysis for Molecular Simulations". arXiv:1909.09830 [physics.comp-ph].

# Compilation

## Dependency Management

### Lua

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

# Nix Usage

We use a deterministic build system to generate both bug reports and uniform
usage statistics.

## Build

Since this project is built with `nix`, we can simply do the following from the
root directory:

```sh
# This will take a long time the first time as it builds the dependencies
nix-build .
# Install into your path
nix-env -if .
# Use anywhere
cd lua_inputs/
yodaStruct -c config.yml
# Use with lua modules
nix-shell --run 'bash' --pure
```

### Language Server Support

To generate a `compile_commands.json` file for working with a language server
like [ccls](https://github.com/MaskRay/ccls) use the following commands:

```sh
# Pure environment
nix-shell --run 'bash' --pure
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=YES ../
cp compile_commands.json ../
```

Note that there is no need to actually compile the project if you simply need to
get the compiler database for the language server.

**Do Not** commit the `.json` file.

### Caveats

Though the build itself is guaranteed to be reproducible as the `nixpkgs` are
also pinned to a particular commit, the `luarocks` dependencies are still local,
since they are determined at runtime.

The above caveats are not relevant when you run it in the shell environment
defined by `shell.nix`. The **recommended usage is to run it in nix-shell**
without **pure**.

```bash
# Install
nix-env -if .
# Go into shell with lfs
nix-shell
# Run the command anywhere
```

## Development

We can simply use the `nix` environment:

```sh
# From the project root
nix-shell
```

# Running

To run the sample inputs, simply move the binary to the project root, or to a
directory where `input/` is a child directory.

```{bash}
# Assuming you are in the build directory
# Check help with -h
./yodaStruct -c ../lua_inputs/config.yml
```

This can also now be tested with a single shell script, which will drop into the
`nix` environment before building and executing the single run test listed
above:

```{bash}
# Just run this
./testBuild.sh
# Or, better yet
nix-build .
nix-env -if .
# If you get a CMake error
nix-build --check .
nix-collect-garbage # then try again [worst case scenario]
```

## Tests

Apart from the [examples](https://docs.dseams.info/pages.html), the test-suite can be run with the `yodaStruct_test` binary.

# Developer Documentation

<!-- TODO: Move this to some other location. -->

## Leaks and performance

While testing for leaks, use `clang` (for
[AddressSanitizer](https://github.com/google/sanitizers/wiki/AddressSanitizer)
and
[LeakSanitizer](https://github.com/google/sanitizers/wiki/AddressSanitizerLeakSanitizer))
and the following:

```{bash}
export CXX=/usr/bin/clang++ && export CC=/usr/bin/clang
cmake .. -DCMAKE_CXX_FLAGS="-pg -fsanitize=address " -DCMAKE_EXE_LINKER_FLAGS=-pg -DCMAKE_SHARED_LINKER_FLAGS=-pg
```

# Contributing

Please ensure that all contributions are formatted according to the
[clang-format](./clang-format) configuration file.

Specifically, consider using the following:

- [Sublime Plugin](https://github.com/rosshemsley/SublimeClangFormat) for users
  of Sublime Text

- [format-all](https://github.com/lassik/emacs-format-all-the-code) for Emacs
- [vim-clang-format](https://github.com/rhysd/vim-clang-format) for Vim
- Visual Studio: http://llvm.org/builds/, or use the [integrated support in Visual Studio 2017](https://blogs.msdn.microsoft.com/vcblog/2018/03/13/clangformat-support-in-visual-studio-2017-15-7-preview-1/)
- Xcode: https://github.com/travisjeffery/ClangFormat-Xcode

Where some of the above suggestions are derived from [this depreciated githook](https://github.com/andrewseidl/githook-clang-format).

Also, do note that we have a `CONTRIBUTING` file you **need to read** to
contribute, for certain reasons, like, common sense.

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
- [Spectra](https://github.com/yixuan/spectra/)
- [Boost Geometry](https://www.boost.org/doc/libs/1_68_0/libs/geometry/doc/html/index.html) for working with different coordinates
- [Boost Math](https://www.boost.org/doc/libs/?view=category_math) for spherical harmonics
- [Blaze](https://bitbucket.org/blaze-lib/blaze/) for very fast modern linear algebra
- [nanoflann](https://github.com/jlblancoc/nanoflann) to calculate nearest neighbors
