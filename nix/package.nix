{
  lib,
  stdenv,
  meson,
  ninja,
  pkg-config,
  eigen,
  blas,
  lapack,
  catch2_3,
  libhwy,
  llvmPackages,
}:

stdenv.mkDerivation (finalAttrs: {
  pname = "seams";
  version = "2.9.2";

  src = lib.fileset.toSource {
    root = ./..;
    fileset = lib.fileset.unions [
      ../meson.build
      ../meson_options.txt
      ../src
      ../tests
      ../input
      ../templates
      ../subprojects/readcon-core.wrap
      ../subprojects/vesin.wrap
      ../subprojects/robin-map.wrap
      ../subprojects/packagefiles
    ];
  };

  nativeBuildInputs = [
    meson
    ninja
    pkg-config
  ];

  buildInputs = [
    eigen
    blas
    lapack
    catch2_3
    libhwy
  ]
  ++ lib.optionals stdenv.cc.isClang [ llvmPackages.openmp ];

  # The meson setup hook already sets -Dwrap_mode=nodownload.
  # auto_features=enabled would force IRA/sphericart/nauty on and fail.
  mesonAutoFeatures = "disabled";
  mesonFlags = [
    (lib.mesonBool "with_tests" true)
    (lib.mesonBool "with_cli" true)
    (lib.mesonBool "with_python" false)
    (lib.mesonEnable "with_lua" false)
    (lib.mesonEnable "with_mpi" false)
    (lib.mesonEnable "with_openmp_offload" false)
  ];

  # Catch2 plus the seams read smoke test. Optional backends (vesin,
  # readcon-core, IRA, sphericart, nauty) stay off unless nixpkgs
  # already provides them; meson skips those probes.
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
