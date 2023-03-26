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

# Define
let
  sources = import ./nix/sources.nix;
  pkgs = import sources.nixpkgs { };
  inherit (pkgs.lib) optional optionals;
  # Import
  buildpkgs = import ./nix { };
  doxygen191 = pkgs.doxygen.overrideAttrs (_: rec {
  name = "doxygen-1.9.1";
  src = pkgs.fetchurl {
    urls = [
      "mirror://sourceforge/doxygen/${name}.src.tar.gz" # faster, with https, etc.
      "http://doxygen.nl/files/${name}.src.tar.gz"
    ];
    sha256 = "1lcif1qi20gf04qyjrx7x367669g17vz2ilgi4cmamp1whdsxbk7";
  };
  });
in pkgs.mkShell {
  # this will make all the build inputs described
  # available to the shell environment
  inputsFrom = [ buildpkgs.yodaStruct pkgs.git ];
  nativeBuildInputs = with pkgs; [ pkgconfig ];
  buildInputs = with pkgs;
    [ gdb ccache ninja git lua luaPackages.luafilesystem doxygen191 which ] ++ optional stdenv.isLinux
    glibcLocales # To allow setting consistent locale on linux
    ++ optional stdenv.isLinux inotify-tools # For file_system
    ++ optional stdenv.isLinux libnotify # For ExUnit
  ;
  # Set up environment vars
  # We unset TERM b/c of https://github.com/NixOS/nix/issues/1056
  shellHook = ''
    stty sane
    export TERM="xterm-256color"
    export LANG="en_US.UTF-8"
    export LC_ALL="en_US.UTF-8"
  '';
}
