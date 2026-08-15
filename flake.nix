{
  description = "d-SEAMS C++ engine (libyodaLib) and the seams CLI";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs =
    { self, nixpkgs }:
    let
      inherit (nixpkgs) lib;
      systems = [
        "x86_64-linux"
        "aarch64-linux"
        "x86_64-darwin"
        "aarch64-darwin"
      ];
      forAllSystems = lib.genAttrs systems;
      pkgsFor = system: nixpkgs.legacyPackages.${system};
    in
    {
      packages = forAllSystems (
        system:
        let
          pkgs = pkgsFor system;
          seams = pkgs.callPackage ./nix/package.nix { };
        in
        {
          inherit seams;
          default = seams;
        }
      );

      apps = forAllSystems (system: {
        default = {
          type = "app";
          program = "${self.packages.${system}.default}/bin/seams";
          meta.description = "Run the seams engine CLI";
        };
      });

      checks = forAllSystems (system: {
        seams = self.packages.${system}.default;
        default = self.checks.${system}.seams;
      });

      devShells = forAllSystems (
        system:
        let
          pkgs = pkgsFor system;
        in
        {
          default = pkgs.mkShell {
            name = "seams-core-dev";
            inputsFrom = [ self.packages.${system}.default ];
            packages = with pkgs; [
              gdb
              ccache
              clang-tools
            ];
          };
        }
      );

      formatter = forAllSystems (system: (pkgsFor system).nixfmt);
    };
}
