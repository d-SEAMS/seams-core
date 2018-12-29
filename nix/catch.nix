# [[file:~/Git/Github/C++/Forks/structureFactor/literateNix.org::*Putting%20Things%20Together][Putting Things Together:3]]
{ stdenv, fetchurl }:

stdenv.mkDerivation rec {
  name = "catch-${version}";
  version = "2.1.2";

  src = fetchurl {
      url = "https://github.com/catchorg/Catch2/releases/download/v2.1.2/catch.hpp";
      sha256 = "e8b8f3109716891aa99b1a8e29cd0d627419bdc4a8d2eeef0d8370aaf8d5e483";
  };

  # It is just the file. No unpacking needed. Seems like we need to create
  # _some_ folder, otherwise we get errors.
  unpackCmd = "mkdir dummy_dir";

  installPhase = ''
    mkdir -p $out/include/catch
    cp ${src} $out/include/catch/catch.hpp
  '';

  meta = {
    description = "A modern, C++-native, header-only, test framework for unit-tests, TDD and BDD - using C++11, C++14, C++17 and later";
    homepage = http://catch-lib.net;
  };
}
# Putting Things Together:3 ends here
