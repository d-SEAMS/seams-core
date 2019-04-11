# [[file:~/Git/Github/C++/Forks/structureFactor/literateNix.org::*Project%20Source][Project Source:1]]
# Usage Example
# nix-shell --argstr compiler gcc5 --run bash
# nix-shell --argstr compiler clang --run bash
# something ? default value ---- Variable declration
# pattern : body ---- Function prototype
{ nixpkgs ? import ./nixpkgs
  , compiler ? "clang"
   }:
  # Define
  let overlay = self: buildpkgs: with buildpkgs; {
  # All the other nix files
  fetchJSON = import ./fetchJSON.nix { inherit (buildpkgs) fetchFromGitHub; };
  catch2 = callPackage ./pkgs/catch2.nix { };
  rang = callPackage ./pkgs/rang.nix { };
  conan = callPackage ./pkgs/conan/conan.nix { };
  fmtlib = callPackage ./pkgs/fmtlib/fmt.nix { };
  yamlCpp = callPackage ./pkgs/yaml-cpp/yaml-cpp.nix { };
  sharkML = callPackage ./pkgs/sharkML/sharkML.nix { };
  # Program expression
    yodaStruct = callPackage ./yodaStruct.nix { };
  }; in
  # Ensure reproducibility
  nixpkgs {
  config = {};
  overlays = [overlay];
  }
# Project Source:1 ends here
