# Note that this script can accept some limited command-line arguments, run
# `julia build_tarballs.jl --help` to see a usage message.
using BinaryBuilder, Pkg

const YGGDRASIL_DIR = "../.."
include(joinpath(YGGDRASIL_DIR, "platforms", "macos_sdks.jl"))

name = "libseams"
version = v"2.14.0"

# seams-core 2.14.0 (release branch head). BinaryBuilder refuses GitHub's
# generated archives, and the project publishes no release tarball, so the
# tree is pinned by commit. Once v2.14.0 is tagged on d-SEAMS/seams-core,
# point the URL there; the commit stays the same.
sources = [
    GitSource("https://github.com/HaoZeke/seams-core.git",
              "18820a224dfebbaf0c327f99f28d8c453f84a221"),
    DirectorySource("./bundled"),
]

# The engine core, its stable C API (libseams), and the C ABI of the R and
# Julia packages (libseams_lang). Every optional backend is off: those come
# from wraps that would be downloaded during the build, or from libraries
# that have no JLL yet.
script = raw"""
cd ${WORKSPACE}/srcdir/seams-core
for p in ${WORKSPACE}/srcdir/patches/*.patch; do
    atomic_patch -p1 "${p}"
done

meson setup build \
    --cross-file="${MESON_TARGET_TOOLCHAIN}" \
    --buildtype=release \
    --wrap-mode=nodownload \
    -Ddefault_library=shared \
    -Dwith_tests=false \
    -Dwith_python=false \
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
    -Dwith_lua=disabled \
    -Dwith_chemfiles=disabled \
    -Dwith_linkcell=disabled \
    -Dwith_minimage=disabled
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
