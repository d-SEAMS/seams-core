# Copyright (c) d-SEAMS developers
# Note that this script can accept some limited command-line arguments, run
# `julia build_tarballs.jl --help` to see a usage message.
using BinaryBuilder, Pkg

const YGGDRASIL_DIR = "../.."
include(joinpath(YGGDRASIL_DIR, "platforms", "macos_sdks.jl"))

name = "libseams"
version = v"2.14.0"

# BinaryBuilder refuses GitHub's generated archives, and seams-core publishes
# no release tarball, so the tree is pinned by commit. v2.14.0 is not tagged.
# v2_14_0_tag_sha is the stand-in for that tag: the release-candidate commit
# on origin/release/2.14-merge. When the tag exists, replace this value with
# `git rev-parse v2.14.0` if the tag is a different commit. The portable
# sources, optional-BLAS meson logic, and non-x86 guards are already in this
# commit, so bundled/patches is not applied.
const v2_14_0_tag_sha = "c1efa8c734214e9117cb77e89258bc9fdd214a6c"

# Pinned to subprojects/backward-cpp.wrap: v1.6's tag commit (the wrap's
# source_hash is a GitHub tarball hash of this same commit; BinaryBuilder
# refuses to take GitHub's generated tag archives directly, see
# check_github_archive in BinaryBuilderBase.jl, so this is a GitSource
# instead), then the same ../../../../subprojects/packagefiles/backward-cpp
# the wrap's own patch_directory applies: meson.build/meson_options.txt/
# LICENSE.build, which upstream backward-cpp does not ship, trimmed of
# backward-cpp's upstream testexe_* targets (see that directory for why).
# with_backward below is on because of these two, not despite nodownload:
# the subproject's meson.build already exists before configure runs, so
# Meson's wrap resolver uses it directly and never needs network.
const backward_cpp_tag_sha = "3bb9240cb15459768adb3e7d963a20e1523a6294"

sources = [
    GitSource("https://github.com/d-SEAMS/seams-core.git", v2_14_0_tag_sha),
    GitSource("https://github.com/bombela/backward-cpp.git", backward_cpp_tag_sha;
               tag = "v1.6", unpack_target = "seams-core/subprojects/backward-cpp-1.6"),
    DirectorySource(joinpath(@__DIR__, "../../../../subprojects/packagefiles/backward-cpp");
                     target = "seams-core/subprojects/backward-cpp-1.6"),
]

# The engine core, its stable C API (libseams), and the C ABI of the R and
# Julia packages (libseams_lang). Every optional backend is off: those come
# from wraps that would be downloaded during the build, or from libraries
# that have no JLL yet. with_backward stays on: the two extra sources
# above pre-seed its subproject directory, so Meson finds a buildable
# backward-cpp already on disk and the nodownload below never applies to it.
script = raw"""
cd ${WORKSPACE}/srcdir/seams-core

meson setup build \
    --cross-file="${MESON_TARGET_TOOLCHAIN}" \
    --buildtype=release \
    --wrap-mode=nodownload \
    -Ddefault_library=shared \
    -Dwith_tests=false \
    -Dwith_cli=true \
    -Dwith_lammps=disabled \
    -Dwith_gromacs=disabled \
    -Dwith_openmm=disabled \
    -Dwith_plumed=disabled \
    -Dwith_mpi=disabled \
    -Dwith_gpulite=disabled \
    -Dwith_openmp_offload=disabled \
    -Dwith_ira=disabled \
    -Dwith_spglib=disabled \
    -Dwith_csm=disabled \
    -Dwith_sphericart=disabled \
    -Dwith_nauty=disabled \
    -Dwith_blas=disabled \
    -Dwith_chemfiles=disabled \
    -Dwith_readcon=disabled \
    -Dwith_linkcell=disabled \
    -Dwith_minimage=disabled \
    -Dwith_spglib=disabled \
    -Dwith_sphericart=disabled \
    -Dwith_nauty=disabled \
    -Dwith_csm=disabled \
    -Dwith_backward=enabled
ninja -C build -j${nproc} src/libyodaLib.${dlext} src/libseams.${dlext} bindings/libseams_lang.${dlext} seams${exeext}

install -Dvm 755 build/src/libyodaLib.${dlext} "${libdir}/libyodaLib.${dlext}"
install -Dvm 755 build/src/libseams.${dlext} "${libdir}/libseams.${dlext}"
install -Dvm 755 build/bindings/libseams_lang.${dlext} "${libdir}/libseams_lang.${dlext}"
install -Dvm 755 build/seams${exeext} "${bindir}/seams${exeext}"
install -Dvm 644 include/seams.h "${includedir}/seams.h"
install -Dvm 644 build/src/seams_version.h "${includedir}/seams_version.h"
install -Dvm 644 bindings/cxx/seams_lang.h "${includedir}/seams/seams_lang.h"
install_license LICENSE
"""

# C++20 headers (<concepts>, <span>, <ranges>) are missing from the default
# macOS SDK.
sources, script = require_macos_sdk("14.5", sources, script; deployment_target = "11.0")

# __int128 in the CRT accumulator rules out 32-bit targets. The Windows and
# FreeBSD ports of the engine are not done yet.
platforms = supported_platforms(; exclude = p -> nbits(p) == 32 || Sys.iswindows(p) || Sys.isfreebsd(p))
platforms = expand_cxxstring_abis(platforms)

products = [
    LibraryProduct("libyodaLib", :libyodaLib),
    LibraryProduct("libseams", :libseams),
    LibraryProduct("libseams_lang", :libseams_lang),
    ExecutableProduct("seams", :seams),
]

dependencies = [
    BuildDependency(PackageSpec(name = "Eigen_jll", uuid = "bc6bbf8a-a594-5541-9c57-10b0d0312c70",
                                version = "3.4")),
    Dependency(PackageSpec(name = "CompilerSupportLibraries_jll", uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae");
               platforms = filter(!Sys.isapple, platforms)),
    Dependency(PackageSpec(name = "LLVMOpenMP_jll", uuid = "1d63c593-3942-5779-bab2-d838dc0a180e");
               platforms = filter(Sys.isapple, platforms)),
]

# C++20 (<numbers>, <span>, std::from_chars on integers) needs GCC 12 and,
# on macOS, LLVM 17.
build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies;
               julia_compat = "1.6", preferred_gcc_version = v"12",
               preferred_llvm_version = v"17")
