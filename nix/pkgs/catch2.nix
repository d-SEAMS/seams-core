# [[file:~/Git/Github/C++/Forks/structureFactor/literateNix.org::*Catch2][Catch2:1]]
{ clangStdenv, fetchurl }:

clangStdenv.mkDerivation rec {
  name = "catch-${version}";
  version = "2.5.0";

  src = fetchurl {
      url = "https://github.com/catchorg/Catch2/releases/download/v2.5.0/catch.hpp";
      sha256 = "a87d5c0417aaf1c3d16565244a1b643e1999d5838d842823731bc18560268f94";
  };

  # This is a header only library. No unpacking needed. Seems like we need to create
  # _some_ folder, otherwise we get errors.
  unpackCmd = "mkdir fake_dir";

  installPhase = ''
    mkdir -p $out/include/catch
    cp ${src} $out/include/catch/catch.hpp
  '';

  meta = {
    description = "A modern, C++-native, header-only, test framework for unit-tests, TDD and BDD - using C++11, C++14, C++17 and later";
    homepage = http://catch-lib.net;
  };
}
# Catch2:1 ends here
