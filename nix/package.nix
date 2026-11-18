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
  llvmPackages,
  git,
  python3,
  cmake,
  util-linux,
  fetchgit,
  fetchurl,
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
}:

let
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
  ovitoWheel = fetchurl {
    url = "https://files.pythonhosted.org/packages/67/99/09355940b7be501f016e481e69c67aa795aac706dd749b0b34ca6f6be482/ovito-3.17.2-cp314-cp314-manylinux_2_34_x86_64.whl";
    sha256 = "1j6zplnpzyn7s87h79cgqkilfpp4ganbb5xv86vwcbz4sdwl8sr3";
  };
  solvisWheel = fetchurl {
    url = "https://files.pythonhosted.org/packages/6c/b4/6664afb3c9aa2747fe2c491adc2b9efbe555d0824f86836711f8b3c9752c/solvis_tools-0.1.5-py3-none-any.whl";
    sha256 = "1s80abs8pz9n9zz9yh04clmyfsnmcxva86q59i6zh2f0rji5cwyx";
  };
  pysideWheel = fetchurl {
    url = "https://files.pythonhosted.org/packages/7b/3c/7597d8dcb36aa8315427635ec6ec9d850ffbad4ce78b9118ac854889f2d3/pyside6_essentials-6.12.0-cp310-abi3-manylinux_2_34_x86_64.whl";
    sha256 = "1n9v7n4nkkjllidcrya0h0ymzzlyv91q4ia8w1d26jdg3pjbbx0n";
  };
  shibokenWheel = fetchurl {
    url = "https://files.pythonhosted.org/packages/49/c2/a04fa003b4360e99d8e3fc2829abacf37396ace63ea1cd616402209ec7c6/shiboken6-6.12.0-cp310-abi3-manylinux_2_34_x86_64.whl";
    sha256 = "0h329n41j6fnzqmmlfg9871aprgdvnnhhf7nqysl8049n1dmdpiv";
  };
  nanobindRoot = "${py.pkgs.nanobind}/lib/python${py.pythonVersion}/site-packages/nanobind";
in

stdenv.mkDerivation (finalAttrs: {
  pname = "seams";
  version = "2.14.0";

  src = lib.fileset.toSource {
    root = ./..;
    fileset = lib.fileset.unions [
      ../meson.build
      ../meson_options.txt
      ../.gitignore
      ../src
      ../include
      ../tests
      ../bindings
      ../engines
      ../bench
      ../gpu_host_kernels
      ../scripts
      ../examples
      ../docs
      ../input
      ../templates
      ../changelog.d
      ../towncrier.toml
      ../README.md
      ../CHANGELOG.md
      ../CITATION.cff
      ../pixi.toml
      ../llms.txt
      ../recipe
      ../data
      ../benchmarks
      ../analysis
      ../Snakefile
      ../cmake
      ../adapters
      ../extras
      ../subprojects/readcon-core.wrap
      ../subprojects/readcon-db.wrap
      ../subprojects/linkcell.wrap
      ../subprojects/minimage.wrap
      ../subprojects/vesin.wrap
      ../subprojects/packagefiles
    ];
  };

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
      freud
      pyvista
      traits
      towncrier
      openmm
    ]))
    openmpi
    lammps
    gromacs
    plumed
  ];

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
    lammps
    gromacs
    openmm
    plumed
  ]
  ++ lib.optionals stdenv.cc.isClang [ llvmPackages.openmp ];

  # The meson setup hook already sets -Dwrap_mode=nodownload.
  # Backends packaged by nixpkgs are required. The others have no package
  # on this pin and cannot use their networked fallbacks in the sandbox.
  mesonFlags = [
    (lib.mesonBool "with_tests" true)
    (lib.mesonBool "with_cli" true)
    (lib.mesonEnable "with_mpi" false)
    (lib.mesonEnable "with_openmp_offload" false)
    (lib.mesonEnable "with_spglib" true)
    (lib.mesonEnable "with_nauty" true)
    (lib.mesonEnable "with_csm" true)
    (lib.mesonEnable "with_lammps" true)
    (lib.mesonEnable "with_gromacs" true)
    (lib.mesonEnable "with_openmm" true)
    (lib.mesonEnable "with_plumed" true)
    (lib.mesonEnable "with_chemfiles" false)
    (lib.mesonEnable "with_readcon" false)
    (lib.mesonEnable "with_ira" false)
    (lib.mesonEnable "with_sphericart" false)
    (lib.mesonEnable "with_gpulite" false)
    (lib.mesonEnable "with_linkcell" false)
    (lib.mesonEnable "with_minimage" false)
  ];

  # Catch2, bindings, command-line examples and available engine backends.
  # docs/check_prose.py and docs/check_tables.py list pages with
  # git ls-files. The filtered source has no repository of its own.
  doCheck = true;
  preCheck = ''
    git -C .. init -q -b main
    git -C .. add -A
    git -C .. -c user.name=d-SEAMS -c user.email=dev@dseams.info commit -q -m sources
    # Scenes invoke "$CC". A bare gcc name resolves to a shim that calls
    # itself. The absolute wrapper does not.
    export CC=$(type -P gcc)
    export CXX=$(type -P g++)
    front=$PWD/pyfront
    rm -rf "$front"
    mkdir -p "$front/dseams" "$front/pydseams"
    cp -a ${pydseamslibSrc}/src/dseams/. "$front/dseams/"
    cp -a ${pydseamslibSrc}/src/pydseams/. "$front/pydseams/"
    chmod -R u+w "$front"
    unzip -q -o ${ovitoWheel} -d "$front"
    unzip -q -o ${solvisWheel} -d "$front"
    unzip -q -o ${pysideWheel} -d "$front"
    unzip -q -o ${shibokenWheel} -d "$front"
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
    export QT_QPA_PLATFORM=offscreen
    "$CXX" -shared -fPIC -std=c++20 -O2 \
      -I${py}/include/python${py.pythonVersion} \
      -I${nanobindRoot}/include \
      -I${robinMap}/include \
      -I../include -Isrc \
      ${pydseamslibSrc}/src/dseams/_core.cpp \
      ${nanobindRoot}/src/nb_combined.cpp \
      -Lsrc -lseams -Wl,-rpath,"$PWD/src" \
      -o "$front/dseams/_core.so"
    "$CXX" -shared -fPIC -std=c++20 -O2 -fopenmp \
      -DSEAMS_HAS_OPENMP=1 -DSEAMS_HAS_HWY=1 \
      -I${py}/include/python${py.pythonVersion} \
      -I${nanobindRoot}/include \
      -I${robinMap}/include \
      -I../include -I../src/include/internal -I../src/include/external \
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
    license = lib.licenses.mit;
    mainProgram = "seams";
    platforms = lib.platforms.linux ++ lib.platforms.darwin;
  };

  passthru = {
    inherit (finalAttrs) version;
  };
})
