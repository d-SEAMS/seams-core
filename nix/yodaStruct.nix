# [[file:~/Git/Github/C++/Forks/structureFactor/literateNix.org::*Expression][Expression:1]]
# Using patterns, and white space negligence
{ clangStdenv
, lua
, conan
, luaPackages
, lib
, boost
, cmake }:
  clangStdenv.mkDerivation {
  name = "yodaStruct";
  src = lib.cleanSource ../.;
  nativeBuildInputs = [
  cmake
  lua
  conan
  ];
  buildInputs = [
  boost
  luaPackages.luafilesystem
  ];
  }
# Expression:1 ends here
