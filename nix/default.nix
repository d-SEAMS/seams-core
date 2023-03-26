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
