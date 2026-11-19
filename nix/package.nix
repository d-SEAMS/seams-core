{
  lib,
  stdenv,
  meson,
  ninja,
  pkg-config,
  fetchurl,
  eigen,
  blas,
  lapack,
  catch2_3,
  libhwy,
  gmp,
  mpfr,
  cxxopts,
  rang,
  llvmPackages,
}:

let
  # Pinned to subprojects/backward-cpp.wrap's source hash. nixpkgs' own
  # backward-cpp installs only a BackwardConfig.cmake meant for
  # add_subdirectory()/add_backward(), not a find_package()-able target
  # (same as conda-forge's, see recipe/build-libseams.sh), so Meson's
  # dependency() cannot locate it and with_backward would otherwise have
  # to fall back to a wrap fetch the sandboxed build phase cannot make.
  # Fetching the upstream tarball here, as a fixed-output derivation,
  # and unpacking it into subprojects/backward-cpp-1.6 before configure,
  # with ../subprojects/packagefiles/backward-cpp (already part of src
  # below) copied on top the same way the wrap's own patch_directory
  # does, reproduces what Meson's wrap resolver would fetch and patch
  # itself, so with_backward stays enabled without any network access
  # during the actual build phase.
  backwardCppSrc = fetchurl {
    url = "https://github.com/bombela/backward-cpp/archive/refs/tags/v1.6.tar.gz";
    sha256 = "c654d0923d43f1cea23d086729673498e4741fb2457e806cfaeaea7b20c97c10";
  };
in
stdenv.mkDerivation (finalAttrs: {
  pname = "seams";
  version = "2.14.0";

  src = lib.fileset.toSource {
    root = ./..;
    fileset = lib.fileset.unions [
      ../meson.build
      ../meson_options.txt
      ../include
      ../cmake
      ../docs
      ../src
      ../engines
      ../bench
      ../gpu_host_kernels
      ../scripts
      ../third_party
      ../tests
      ../bindings
      ../benchmarks
      ../input
      ../templates
      ../subprojects/readcon-core.wrap
      ../subprojects/readcon-db.wrap
      ../subprojects/linkcell.wrap
      ../subprojects/minimage.wrap
      ../subprojects/vesin.wrap
      ../subprojects/cxxopts.wrap
      ../subprojects/backward-cpp.wrap
      ../subprojects/rang.wrap
      ../subprojects/packagefiles
    ];
  };

  nativeBuildInputs = [
    meson
    ninja
    pkg-config
  ];

  # Reproduces what Meson's wrap resolver does for backward-cpp.wrap: the
  # upstream v1.6 tarball, then subprojects/packagefiles/backward-cpp's
  # meson.build/meson_options.txt/LICENSE.build overlaid on top (upstream
  # backward-cpp ships none of those three; see that directory for why
  # its meson.build is trimmed of backward-cpp's own testexe_* targets).
  # Doing it here means the subproject's meson.build already exists
  # before configure runs, so Meson's wrap resolver uses it directly and
  # never needs network.
  postPatch = ''
    mkdir -p subprojects/backward-cpp-1.6
    tar -xzf ${backwardCppSrc} --strip-components=1 -C subprojects/backward-cpp-1.6
    cp subprojects/packagefiles/backward-cpp/* subprojects/backward-cpp-1.6/
  '';

  buildInputs = [
    eigen
    blas
    lapack
    catch2_3
    libhwy
    gmp
    mpfr
    cxxopts
    rang
  ]
  ++ lib.optionals stdenv.cc.isClang [ llvmPackages.openmp ];

  # The meson setup hook already sets -Dwrap_mode=nodownload.
  # auto_features=enabled would force IRA/sphericart/nauty on and fail;
  # with_backward is explicitly re-enabled below despite that, since
  # postPatch has already made its subproject buildable offline.
  mesonAutoFeatures = "disabled";
  mesonFlags = [
    (lib.mesonBool "with_tests" true)
    (lib.mesonBool "with_cli" true)
    (lib.mesonEnable "with_mpi" false)
    (lib.mesonEnable "with_openmp_offload" false)
    (lib.mesonEnable "with_backward" true)
  ];

  # Catch2 plus the seams read smoke test. Optional backends (vesin,
  # readcon-core, IRA, sphericart, nauty) stay off unless nixpkgs already
  # provides them; meson skips those probes. backward is on: its
  # subproject is pre-seeded by postPatch above.
  doCheck = true;

  meta = {
    description = "d-SEAMS C++ engine and seams CLI";
    homepage = "https://dseams.info";
    license = lib.licenses.mit;
    mainProgram = "seams";
    platforms = lib.platforms.linux ++ lib.platforms.darwin;
  };

  passthru = {
    inherit (finalAttrs) version;
  };
})
