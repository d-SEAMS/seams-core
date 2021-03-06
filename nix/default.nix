#-----------------------------------------------------------------------------------
# d-SEAMS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# A copy of the GNU General Public License is available at
# http://www.gnu.org/licenses/
#-----------------------------------------------------------------------------------

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
  # Package for testing
  catch2 = callPackage ./pkgs/catch2.nix { };
  # Package for testing
  fmtlib = callPackage ./pkgs/fmtlib/fmt.nix { };
  # Package for testing
  yamlCpp = callPackage ./pkgs/yaml-cpp/yaml-cpp.nix { };
  # Program expression
    yodaStruct = callPackage ./yodaStruct.nix { };
  }; in
  # Ensure reproducibility
  nixpkgs {
  config = {};
  overlays = [overlay];
  }
# Project Source:1 ends here
