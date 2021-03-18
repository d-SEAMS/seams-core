# Define
let
  sources = import ./nix/sources.nix;
  pkgs = import sources.nixpkgs { };
  inherit (pkgs.lib) optional optionals;
  # Import
  buildpkgs = import ./nix { };
  doxygen = pkgs.doxygen.overrideAttrs (_: { version = "1.9.1"; });
in pkgs.mkShell {
  # this will make all the build inputs described
  # available to the shell environment
  inputsFrom = [ buildpkgs.yodaStruct pkgs.git ];
  nativeBuildInputs = with pkgs; [ pkgconfig ];
  buildInputs = with pkgs;
    [ gdb git lua luaPackages.luafilesystem doxygen ] ++ optional stdenv.isLinux
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
