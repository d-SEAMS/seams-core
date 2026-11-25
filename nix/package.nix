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
  gmp,
  mpfr,
  cxxopts,
  rang,
  llvmPackages,
  git,
  python3,
  cmake,
  util-linux,
  fetchgit,
  fetchurl,
  fetchzip,
  rustPlatform,
  rustc,
  cargo,
  robinMap,
  unzip,
  libglvnd,
  glib,
  brotli,
  libxkbcommon,
  krb5,
  xorg,
  fontconfig,
  vulkanLoader,
  dbus,
  spglib,
  nauty,
  openbabel,
  lammps,
  gromacs,
  openmm,
  plumed,
  openmpi,
  fftwFloat,
  runCommand,
  gnutar,
  gzip,
}:

let
  available = lib.meta.availableOn stdenv.hostPlatform;
  withLammps = available lammps;
  withGromacs = available gromacs;
  withOpenmm = available openmm;
  withPlumed = available plumed;
  lammpsWithPlugin = lammps.overrideAttrs (old: {
    cmakeFlags = (old.cmakeFlags or [ ]) ++ [ "-DPKG_PLUGIN=ON" ];
  });
  engineInputs = lib.optional withLammps lammpsWithPlugin
    ++ lib.optional withGromacs gromacs
    ++ lib.optional withPlumed plumed;
  testWheels = import ./test-wheels.nix {
    inherit fetchurl;
    system = stdenv.hostPlatform.system;
  };
  bindingSource = name: rev: hash: {
    inherit rev;
    src = fetchzip {
      url = "https://codeload.github.com/d-SEAMS/${name}/tar.gz/${rev}";
      extension = "tar.gz";
      inherit hash;
    };
  };
  bindingSources = {
    PydSEAMSlib = bindingSource "PydSEAMSlib" "115f0ca9df37107c8ca824088b5d813fe41b2c3e"
      "sha256-+flLctn4JkA/qNS999KGCf6SIYJ9vcEPdZcxNnKBJMY=";
    "Dseams.jl" = bindingSource "Dseams.jl" "24bd4e4d28956fc69a0bf69a86e36145505e7ffe"
      "sha256-ZbwL/Qd00IZ1PJ9PrC9vkwRzjg2PaSXPamd5TfbqgtY=";
    dseams-r = bindingSource "dseams-r" "25dfa13e1b1542d88ac768ecd053ebcfe10c72da"
      "sha256-3KF4IwNGnHCEfGOfCvyMTj2KvdRYk4k4V4ITZ1b7kRk=";
    dseams-rs = bindingSource "dseams-rs" "c462541512940a3101485a0f731cd1b3e6883fcd"
      "sha256-Iedu4qlFRDMX0roxTwBxoaX5cLotkIpkaVThLidee+o=";
    dseams-fortran = bindingSource "dseams-fortran" "aafd4b85440c616338e2831ba3bdcdf64e16b2d2"
      "sha256-RwV+wyDaKDcQYIXeG1YRUNomcU8LrjBFgvYNvvIm8v0=";
    yodaStruct = bindingSource "yodaStruct" "79b6daa1ec74f37451765f537a0d0c37407380fe"
      "sha256-TlC4TAI3LQB95G1+oZRdChp7L+hHWjRWi9mBonnqdg0=";
  };
  # freud needs parsnip. The wheel's distribution name is parsnip_cif.
  py = python3.override {
    packageOverrides = self: super: {
      parsnip = super.parsnip.overridePythonAttrs (_: {
        # The wheel is named parsnip_cif, so the metadata hook looks up
        # the wrong distribution. The import check still runs.
        dontCheckPythonMetadata = true;
      });
      freud = super.freud.overridePythonAttrs (_: {
        dontCheckPythonMetadata = true;
      });
    };
  };
  pydseamslibSrc = fetchgit {
    url = "https://github.com/d-SEAMS/PydSEAMSlib.git";
    rev = "e2bdf66e671b0cc527b5939cf1b3acfd01ebe727";
    hash = "sha256-Too0mXaa+ISNcMA0k/+MdKwHE3i2Lhj5OrMjL4mCoUI=";
  };
  solvisWheel = fetchurl {
    url = "https://files.pythonhosted.org/packages/6c/b4/6664afb3c9aa2747fe2c491adc2b9efbe555d0824f86836711f8b3c9752c/solvis_tools-0.1.5-py3-none-any.whl";
    sha256 = "1s80abs8pz9n9zz9yh04clmyfsnmcxva86q59i6zh2f0rji5cwyx";
  };
  # This pin ships towncrier 25, which joins a Markdown title to the notes
  # with no blank line. The changelog was rendered by towncrier 23.11, the
  # release the pixi lock still runs, so the check uses that release.
  towncrierWheel = fetchurl {
    url = "https://files.pythonhosted.org/packages/36/5a/6546c8ea2bee64c9e9552e06fdb31118a4f23a12aa161581b5fefd070152/towncrier-23.11.0-py3-none-any.whl";
    hash = "sha256-LlGcphlCbRiePJjJlVj+i+UMnO0T6h/CCko1OpXS3tc=";
  };
  nucleationTrajectory = fetchurl {
    url = "https://ndownloader.figshare.com/files/20418309";
    hash = "sha256-w8dQULt3Ys8NNBfi+6syZCLyPqBKZLb2nq8zM05NFAM=";
  };
  brineFrameEarly = fetchurl {
    url = "https://huggingface.co/datasets/d-SEAMS/brine-front-230K/resolve/4828e33ecf419c70643b2d466bfa4f283b4138dc/data/frames/T230_m20_r0_t0.xyz";
    hash = "sha256-3FWDL8rmZKpeUoKDs76Idi798dq/TbFSC7HIBkAbYSg=";
  };
  brineFrameLate = fetchurl {
    url = "https://huggingface.co/datasets/d-SEAMS/brine-front-230K/resolve/4828e33ecf419c70643b2d466bfa4f283b4138dc/data/frames/T230_m20_r0_tend.xyz";
    hash = "sha256-3LwF+B8yG0UZXSWc1MIuaTQPnBrlfZPEenTJlcWN0F0=";
  };
  gromacsLegacyHeaders = runCommand "gromacs-${gromacs.version}-source" {
    nativeBuildInputs = [
      gnutar
      gzip
    ];
  } ''
    mkdir -p "$out"
    tar -xzf ${gromacs.src} --strip-components=1 -C "$out"
  '';
  minimageSrc = fetchzip {
    url = "https://github.com/lode-org/minimage/archive/refs/tags/v0.1.4.tar.gz";
    hash = "sha256-8Y9UFgj7hEiXcRNpUa0cHx6aT01QomJWCvvAuBfXya0=";
  };
  linkcellSrc = fetchzip {
    url = "https://github.com/d-SEAMS/linkcell/archive/refs/tags/v0.6.0.tar.gz";
    hash = "sha256-C1RPSOQHsBjfYxsw5Si6Q4+B+YGlQjNHitinndAv8XA=";
  };
  vesinSrc = fetchzip {
    url = "https://github.com/Luthaf/vesin/archive/refs/tags/v0.5.7.tar.gz";
    hash = "sha256-JsU+tDyATW6wmix7Tjwa5RnLRjB6b+sYSkMJIt5M5Yo=";
  };
  linkcellCargoVendor = rustPlatform.fetchCargoVendor {
    src = linkcellSrc;
    hash = "sha256-2QscNVUPmCrwO1cAHWnDZAtU8Hsr6Bl/WTO65WwDsC8=";
  };
  minimageCargoVendor = rustPlatform.fetchCargoVendor {
    src = minimageSrc;
    hash = "sha256-8NSCv17ACkQbkBnycIZ3K9mUUOWwYX4u7OMii2F2znk=";
  };
  cargoHome = runCommand "seams-cargo-home" { } ''
    mkdir -p "$out/source-registry-0"
    cp -R ${linkcellCargoVendor}/source-registry-0/. "$out/source-registry-0/"
    chmod -R u+w "$out/source-registry-0"
    cp -R --no-clobber ${minimageCargoVendor}/source-registry-0/. "$out/source-registry-0/"
    cat > "$out/config.toml" <<EOF
[source.vendored-source-registry-0]
directory = "$out/source-registry-0"

[source.crates-io]
replace-with = "vendored-source-registry-0"
EOF
  '';
  nanobindRoot = "${py.pkgs.nanobind}/lib/python${py.pythonVersion}/site-packages/nanobind";
  backwardCppSrc = fetchurl {
    url = "https://github.com/bombela/backward-cpp/archive/refs/tags/v1.6.tar.gz";
    sha256 = "c654d0923d43f1cea23d086729673498e4741fb2457e806cfaeaea7b20c97c10";
  };
in
stdenv.mkDerivation (finalAttrs: {
  pname = "seams";
  version =
    let
      afterProject = builtins.elemAt (lib.splitString "\nproject(" ("\n" + builtins.readFile ../meson.build)) 1;
      projectBlock = builtins.head (lib.splitString ")" afterProject);
      lines = lib.splitString "\n" projectBlock;
      matchLine =
        line:
        builtins.match "[ \t]*version:[ \t]*'([0-9]+\\.[0-9]+\\.[0-9]+)',?[ \t]*" line;
      matches = lib.filter (m: m != null) (map matchLine lines);
    in
    if matches == [ ] then
      throw "Could not extract project version from project() block in ../meson.build"
    else
      builtins.head (builtins.head matches);

  src = lib.cleanSource ../.;

  postPatch = ''
    mkdir -p subprojects/backward-cpp-1.6
    tar -xzf ${backwardCppSrc} --strip-components=1 -C subprojects/backward-cpp-1.6
    cp subprojects/packagefiles/backward-cpp/* subprojects/backward-cpp-1.6/
    patchShebangs --build scripts
    mkdir -p docs/.cache/bindings
    ${lib.concatStringsSep "\n" (lib.mapAttrsToList (name: source: ''
      mkdir -p docs/.cache/bindings/${name}
      cp -R ${source.src}/. docs/.cache/bindings/${name}/
      printf '%s\n' '${source.rev}' > docs/.cache/bindings/${name}/.seams-vendored-rev
    '') bindingSources)}
    cp -R ${minimageSrc} subprojects/minimage
    cp -R ${linkcellSrc} subprojects/linkcell
    cp -R ${vesinSrc} subprojects/vesin
    chmod -R u+w subprojects/minimage subprojects/linkcell subprojects/vesin
    patch -d subprojects/linkcell -p1 < subprojects/packagefiles/linkcell-omp-test.patch
    cp -R subprojects/packagefiles/vesin/. subprojects/vesin/
    patch -d subprojects/vesin -p1 < subprojects/packagefiles/vesin-null-reset.patch
    mkdir -p docs/.cache/figshare-11448702
    cp ${nucleationTrajectory} docs/.cache/figshare-11448702/nucleation.lammpstrj
    mkdir -p docs/.cache/brine-front-230K/4828e33ecf419c70643b2d466bfa4f283b4138dc
    cp ${brineFrameEarly} docs/.cache/brine-front-230K/4828e33ecf419c70643b2d466bfa4f283b4138dc/T230_m20_r0_t0.xyz
    cp ${brineFrameLate} docs/.cache/brine-front-230K/4828e33ecf419c70643b2d466bfa4f283b4138dc/T230_m20_r0_tend.xyz
  '';

  nativeBuildInputs = [
    meson
    ninja
    pkg-config
    git
    cmake
    util-linux
    # Meson treats the interpreter it is running as python3. A Python
    # on PATH is the one the test suite actually executes.
    unzip
    (py.withPackages (ps: with ps; [
      numpy
      pillow
      matplotlib
      scipy
      ase
      pytest
      freud
      pyvista
      traits
      towncrier
    ] ++ lib.optional withOpenmm ps.openmm))
    rustc
    cargo
  ] ++ lib.optional (available openmpi) openmpi ++ engineInputs;


  buildInputs = [
    eigen
    blas
    lapack
    catch2_3
    libhwy
    gmp
    mpfr
    spglib
    nauty
    openbabel
    fftwFloat
    cxxopts
    rang
  ]
  ++ engineInputs
  ++ lib.optional withOpenmm openmm
  ++ lib.optional stdenv.cc.isClang llvmPackages.openmp;

  # The meson setup hook already sets -Dwrap_mode=nodownload.
  # Engine integrations use the packages available for the host platform.
  # The sandbox cannot use networked fallbacks for missing dependencies.
  mesonFlags = [
    (lib.mesonBool "with_tests" true)
    (lib.mesonBool "with_cli" true)
    (lib.mesonEnable "with_mpi" false)
    (lib.mesonEnable "with_openmp_offload" false)
    (lib.mesonEnable "with_backward" true)
    (lib.mesonEnable "with_spglib" true)
    (lib.mesonEnable "with_nauty" true)
    (lib.mesonEnable "with_csm" true)
    (lib.mesonEnable "with_lammps" withLammps)
    (lib.mesonEnable "with_gromacs" withGromacs)
    (lib.mesonEnable "with_openmm" withOpenmm)
    (lib.mesonEnable "with_plumed" withPlumed)
    (lib.mesonEnable "with_chemfiles" false)
    (lib.mesonEnable "with_readcon" false)
    (lib.mesonEnable "with_ira" false)
    (lib.mesonEnable "with_sphericart" false)
    (lib.mesonEnable "with_gpulite" false)
    (lib.mesonEnable "with_linkcell" true)
    (lib.mesonEnable "with_minimage" true)
  ]
  ++ lib.optional withLammps (lib.mesonOption "lammps_src_include" "${lammps.src}/src")
  ++ lib.optional withGromacs (lib.mesonOption "gromacs_legacy_include" "${gromacsLegacyHeaders}/api/legacy/include");

  preConfigure = ''
    export CARGO_HOME=${cargoHome}
    export CARGO_NET_OFFLINE=true
  '';

  # Catch2, bindings, command-line examples and available engine backends.
  # docs/check_prose.py and docs/check_tables.py list pages with
  # git ls-files. The filtered source has no repository of its own.
  doCheck = true;
  preCheck = ''
    export SEAMS_NO_NETWORK=1
    git -C .. init -q -b main
    git -C .. add -A
    git -C .. -c user.name=d-SEAMS -c user.email=dev@dseams.info commit -q -m sources
    # Scene tests need the absolute compiler wrappers for this platform.
    export CC=${stdenv.cc}/bin/cc
    export CXX=${stdenv.cc}/bin/c++
    ${lib.optionalString withOpenmm "export OPENMM_DIR=${openmm}"}
    front=$PWD/pyfront
    rm -rf "$front"
    mkdir -p "$front/dseams" "$front/pydseams"
    cp -a ${pydseamslibSrc}/src/dseams/. "$front/dseams/"
    cp -a ${pydseamslibSrc}/src/pydseams/. "$front/pydseams/"
    chmod -R u+w "$front"
    unzip -q -o ${testWheels.ovito} -d "$front"
    unzip -q -o ${solvisWheel} -d "$front"
    unzip -q -o ${testWheels.pyside} -d "$front"
    unzip -q -o ${testWheels.shiboken} -d "$front"
    unzip -q -o ${towncrierWheel} -d "$front"
    ${lib.optionalString stdenv.hostPlatform.isLinux ''
      export LD_LIBRARY_PATH=${lib.makeLibraryPath [
      libglvnd
      glib
      brotli
      libxkbcommon
      krb5
      stdenv.cc.cc.lib
      xorg.libX11
      fontconfig
      vulkanLoader
      dbus
      ]}:$LD_LIBRARY_PATH
    ''}
    export QT_QPA_PLATFORM=offscreen
    "$CXX" -shared -fPIC -std=c++20 -O2 ${lib.optionalString stdenv.hostPlatform.isDarwin "-Wl,-undefined,dynamic_lookup"} \
      -I${py}/include/python${py.pythonVersion} \
      -I${nanobindRoot}/include \
      -I${robinMap}/include \
      -I../include -Isrc \
      ${pydseamslibSrc}/src/dseams/_core.cpp \
      ${nanobindRoot}/src/nb_combined.cpp \
      -Lsrc -lseams -Wl,-rpath,"$PWD/src" \
      -o "$front/dseams/_core.so"
    "$CXX" -shared -fPIC -std=c++20 -O2 -fopenmp ${lib.optionalString stdenv.hostPlatform.isDarwin "-Wl,-undefined,dynamic_lookup"} \
      -DSEAMS_HAS_OPENMP=1 -DSEAMS_HAS_HWY=1 \
      -I${py}/include/python${py.pythonVersion} \
      -I${nanobindRoot}/include \
      -I${robinMap}/include \
      -I../include -I../src/include/internal -I../third_party/argum \
      -I../subprojects/linkcell/include \
      -I${libhwy}/include \
      -isystem ${eigen}/include/eigen3 \
      -Isrc \
      ${pydseamslibSrc}/src/bindings.cpp \
      ${nanobindRoot}/src/nb_combined.cpp \
      -Lsrc -lyodaLib -lseams -Wl,-rpath,"$PWD/src" \
      -o "$front/pydseams/yoda.so"
    export PYTHONPATH="$front''${PYTHONPATH:+:}$PYTHONPATH"
  '';

  meta = {
    description = "d-SEAMS C++ engine and seams CLI";
    homepage = "https://dseams.info";
    # The CLI links third_party/argum (BSD-3-Clause). csm is off here.
    license = with lib.licenses; [ mit bsd3 ];
    mainProgram = "seams";
    platforms = lib.platforms.linux ++ lib.platforms.darwin;
  };

  passthru = {
    inherit (finalAttrs) version;
  };
})
