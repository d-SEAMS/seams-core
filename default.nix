# [[file:~/Git/Github/C++/Forks/structureFactor/literateNix.org::finalNix][finalNix]]
# Nix skeleton for compiler, cmake, boost.
# Dependencies (boost and the rest) are built with selected compiler (for ABI compatibility).
# Examples:
#   nix-shell --argstr compiler gcc5 --run 'mkdir build && cd build && cmake .. && cmake --build .'
#   nix-shell --argstr compiler gcc6 --run 'mkdir build && cd build && cmake .. && cmake --build .'
#   nix-shell --argstr compiler clang_38 --run 'mkdir build && cd build && cmake .. && cmake --build .'
{ nixpkgs ? import <nixpkgs> {}, compiler ? "clang" }:
let
  stdenv = nixpkgs.overrideCC nixpkgs.stdenv nixpkgs.${compiler};
  in rec {
    myProject = stdenv.mkDerivation {
      name = "supaaYoda";
      version = "dev-0.1";
      nativeBuildInputs = [ nixpkgs.boost nixpkgs.lua nixpkgs.cmake nixpkgs.luarocks ];
      buildInputs = [  ];
  };
  }
# finalNix ends here
