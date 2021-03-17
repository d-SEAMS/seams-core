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
, gsl
, cmake }:
  clangStdenv.mkDerivation {
  name = "yodaStruct";
  src = lib.cleanSource ../.;
  enableParallelBuilding=true;
  nativeBuildInputs = [
  fmtlib
  cmake
  lua
  luaPackages.luafilesystem
  ];
  buildInputs = [
  yamlCpp
  eigen
  catch2
  boost
  gsl
  liblapack
  blas
  lua
  luaPackages.luafilesystem
  ];
  }
# Expression:1 ends here
