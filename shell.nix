{ pkgs ? import <nixpkgs> {} }:
  # Define
  let
    inherit (pkgs.lib) optional optionals;
  # Import
    buildpkgs = import ./nix {};
  in
pkgs.mkShell {
  # this will make all the build inputs from hello and gnutar
  # available to the shell environment
  inputsFrom =  [
    buildpkgs.yodaStruct
    pkgs.git
                ];
  buildInputs = with pkgs; [
    git
  ]
  ++ optional stdenv.isLinux glibcLocales # To allow setting consistent locale on linux
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
