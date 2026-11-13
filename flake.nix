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

      apps = forAllSystems (
        system:
        let
          pkgs = pkgsFor system;
          # One app per paper table. Run from the checkout: nix run .#<name>
          tables = {
            sitting = "benchmarks/paper-cells/sitting.csv";
            "tab-rings-214" = "benchmarks/paper-cells/tab-rings-214.csv";
            "tab-ql-eq" = "benchmarks/paper-cells/tab-ql-eq.csv";
            "tab-neighbours" = "benchmarks/paper-cells/tab-neighbours.csv";
            "tab-chill" = "benchmarks/paper-cells/tab-chill.csv";
            "tab-cages" = "benchmarks/paper-cells/tab-cages.csv";
            "tab-cli" = "benchmarks/paper-cells/tab-cli.csv";
            "tab-lit" = "benchmarks/paper-cells/tab-lit.csv";
            "tab-hop-stride" = "benchmarks/paper-cells/tab-hop-stride.csv";
            "tab-sym-eq" = "benchmarks/paper-cells/tab-sym-eq.csv";
            "symmetry-vs-ptm" = "benchmarks/paper-cells/symmetry-vs-ptm.csv";
            "ptm-ice" = "benchmarks/paper-cells/ptm-ice.csv";
          };
          one =
            name: target: {
              type = "app";
              program = "${pkgs.writeShellScript "paper-${name}" ''
                set -eu
                if [ ! -f Snakefile ]; then
                  echo "run this from the seams-core checkout" >&2
                  exit 1
                fi
                exec ${pkgs.snakemake}/bin/snakemake -j2 ${target}
              ''}";
              meta.description = "Regenerate ${target}";
            };
        in
        {
          default = {
            type = "app";
            program = "${self.packages.${system}.default}/bin/seams";
            meta.description = "Run the seams engine CLI";
          };
        }
        // lib.mapAttrs one tables
      );

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
              snakemake
            ];
          };
        }
      );

      formatter = forAllSystems (system: (pkgsFor system).nixfmt);
    };
}
