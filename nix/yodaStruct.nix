#-----------------------------------------------------------------------------------
# d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
#
# Copyright (c) 2018--present d-SEAMS core team
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the MIT License as published by
# the Open Source Initiative.
#
# A copy of the MIT License is included in the LICENSE file of this repository.
# You should have received a copy of the MIT License along with this program.
# If not, see <https://opensource.org/licenses/MIT>.
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
  liblapack
  blas
  lua
  luaPackages.luafilesystem
  ];
  }
# Expression:1 ends here
