#-----------------------------------------------------------------------------------
# d-SEAMS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# A copy of the GNU General Public License is available at
# http://www.gnu.org/licenses/
#-----------------------------------------------------------------------------------

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
