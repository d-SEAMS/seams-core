# [[file:~/Git/Github/C++/Forks/structureFactor/literateNix.org::*Expression][Expression:1]]
# Using patterns, and white space negligence
{ clangStdenv
, catch2
, fmtlib
, yamlCpp
, sharkML
, lua
, luaPackages
, lib
, rang
, boost
, cmake
, prod ? true}:
  clangStdenv.mkDerivation {
  name = "yodaStruct";
  src = lib.cleanSource ../.;
  nativeBuildInputs = [
  catch2
  fmtlib
  cmake
  rang
  lua
  ];
  buildInputs = [
  yamlCpp
  sharkML
  boost
  luaPackages.luafilesystem
  ];
  cmakeFlags = [
    (lib.optional prod "-DCMAKE_EXPORT_COMPILE_COMMANDS=OFF")
];
    # postInstall =  o.postInstall + ''
    #   source $stdenv/setup
    #   cp yodaStruct $out/bin
    #   cp -r lib $out
    # '';
  }
# Expression:1 ends here
