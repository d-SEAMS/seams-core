# d-SEAMS

**Deferred Structural Elucidation Analysis for Molecular Simulations**

[![Build Status](https://travis-ci.org/d-SEAMS/seams-core.svg?branch=master)](https://travis-ci.org/d-SEAMS/seams-core)
[![built with nix](https://builtwithnix.org/badge.svg)](https://builtwithnix.org)

Check our build status [here](https://travis-ci.org/d-SEAMS/seams-core).
The docs themselves are [here](https://docs.dseams.info) and
development is ongoing [on GitHub](https://github.com/d-SEAMS/seams-core). We
also have [a Zenodo community](https://zenodo.org/communities/d-seams/) for user-contributions like reviews, testimonials
and tutorials. Trajectories are hosted [on figshare](https://figshare.com/projects/d-SEAMS_Datasets/73545).

\brief The C++ core of d-SEAMS, a molecular dynamics trajectory analysis engine.

\note The <a href="pages.html">wiki</a> describes the examples and how to obtain
the data-sets (trajectories) <a
href="https://figshare.com/projects/d-SEAMS_Datasets/73545">from figshare</a>.

\warning **If** you are unwilling to use the `nix` build system, then **please note** that you must manage the dependencies MANUALLY, including the compiler versions.

# Citation

This software is being actively developed and written up. If you do use it in an
academic capacity, for now please cite [the following preprint](https://arxiv.org/abs/1909.09830):

    Goswami, R.; Goswami, A.; Singh, J. K. (2019). "d-SEAMS: Deferred Structural Elucidation Analysis for Molecular Simulations". arXiv:1909.09830 [physics.comp-ph].

# Compilation with Nix

We use a deterministic build system to generate both bug reports and uniform
usage statistics. This also handles the `lua` scripting engine.

\note The lua functions are documented on the [wiki](https://docs.dseams.info/md_markdown_luafunctions)

## Build

Since this project is built with `nix`, we can simply do the following from the
root directory (longer method):

```sh
# Make sure there are no artifacts
rm -rf build
# This will take a long time the first time as it builds the dependencies
nix-build . # Optional
# Install into your path
nix-env -if . # Required
# Run the command anywhere
yodaStruct -c lua_inputs/config.yml
```

A faster method of building the software is by using the [cachix binary cache](https://dseams.cachix.org/) as shown:

```bash
# Install cachix
nix-env -iA cachix -f https://cachix.org/api/v1/install
# Use the binary cache
cachix use dseams
# Faster with the cache than building from scratch
nix-build . # Optional
# Install into your path
nix-env -if . # Required
# Run the command anywhere
yodaStruct -c lua_inputs/config.yml
```

\note The paths in the `.yml` should be **relative to the folder from which the binary is called**.

If you're confused about how to handle the relative paths, run the command `yodaStruct -c lua_inputs/config.yml` in the top-level directory, and set the paths relative to the top-level directory. This is the convention used in the examples as well.

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

## Development

We can simply use the `nix` environment:

```sh
# From the project root
nix-shell
```

# Running

This is built completely with nix:

```{bash}
# Install systemwide
nix-env -if .
```

To run the sample inputs, simply install the software, and ensure that `input/` is a child directory.

```{bash}
# Assuming you are in the src directory
# Check help with -h
yodaStruct -c lua_inputs/config.yml
```

## Tests

Apart from the [examples](https://docs.dseams.info/pages.html), the test-suite
can be run with the `yodaStruct_test` binary, which will drop into the
`nix` environment before building and executing `gdb`:

```{bash}
# Just run this
./testBuild.sh
# quit gdb with quit
# Go run the test binary
cd shellBuild
./yodaStruct_test
```

Do note that the regular installation via `nix-env` runs the tests before the installation

# Developer Documentation

<!-- TODO: Move this to some other location. -->

Test the build with nix:

```bash
nix-build .
# Outputs are in ./result
# If you get a CMake error
rm -rf build
nix-store --delete /nix/store/$whatever # $whatever is the derivation complaining
nix-collect-garbage # then try again [worst case scenario]
```

## Leaks and performance

While testing for leaks, use `clang` (for
[AddressSanitizer](https://github.com/google/sanitizers/wiki/AddressSanitizer)
and
[LeakSanitizer](https://github.com/google/sanitizers/wiki/AddressSanitizerLeakSanitizer))
and the following:

```{bash}
# From the developer shell
export CXX=/usr/bin/clang++ && export CC=/usr/bin/clang
cmake .. -DCMAKE_CXX_FLAGS="-pg -fsanitize=address " -DCMAKE_EXE_LINKER_FLAGS=-pg -DCMAKE_SHARED_LINKER_FLAGS=-pg
```

# Overview

As of Mon Jan 20 15:57:18 2020, the lines of code calculated by
[cloc](http://cloc.sourceforge.net/) are as follows:

![Cloc Lines](images/cloc-2020-01-20_15-56.png)

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
