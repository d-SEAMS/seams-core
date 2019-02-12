# [[file:~/Git/Github/C++/Forks/structureFactor/literateNix.org::*Expression][Expression:1]]
# Using patterns, and white space negligence
{ clangStdenv
, catch2
, fmtlib
, yamlCpp
, eigen
, lua
, luaPackages
, liblapack
, blas
, lib
, boost
, cmake }:
  clangStdenv.mkDerivation {
  name = "yodaStruct";
  src = lib.cleanSource ../.;
  nativeBuildInputs = [
  fmtlib
  cmake
  lua
  ];
  buildInputs = [
  yamlCpp
  eigen
  catch2
  boost
  liblapack
  blas
  luaPackages.luafilesystem
  ];
  }
# Expression:1 ends here
