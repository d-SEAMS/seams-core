#-----------------------------------------------------------------------------------
# d-SEAMS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# A copy of the GNU General Public License is available at
# http://www.gnu.org/licenses/
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
    [ gdb vim ccache ninja git lua luaPackages.luafilesystem doxygen191 which ] ++ optional stdenv.isLinux
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
