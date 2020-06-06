# [[file:~/Git/Github/C++/Forks/structureFactor/literateNix.org::*Catch2][Catch2:1]]
{ clangStdenv, fetchurl }:

clangStdenv.mkDerivation rec {
  name = "catch-${version}";
  version = "2.6.0";

  src = fetchurl {
      url = "https://github.com/catchorg/Catch2/releases/download/v2.6.0/catch.hpp";
      sha256 = "1zkm95z5f4ih0fs0fna7ya1ryy39lhg9hwbirzhvc8a79nrk6qd8";
  };

  # This is a header only library. No unpacking needed. Seems like we need to create
  # _some_ folder, otherwise we get errors.
  unpackCmd = "mkdir fake_dir";

  installPhase = ''
    mkdir -p $out/include/catch2
    cp ${src} $out/include/catch2/catch.hpp
  '';

  meta = {
    description = "A modern, C++-native, header-only, test framework for unit-tests, TDD and BDD - using C++11, C++14, C++17 and later";
    homepage = http://catch-lib.net;
  };
}
# Catch2:1 ends here
