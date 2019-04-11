{ pkgs ? import <nixpkgs> {} }:
  # Define
  let
  # Import
    buildpkgs = import ./nix {};
  in
pkgs.mkShell {
  # this will make all the build inputs from hello and gnutar
  # available to the shell environment
  inputsFrom =  [ buildpkgs.yodaStruct ];
  # buildInputs = [ pkgs.hello ];
}
