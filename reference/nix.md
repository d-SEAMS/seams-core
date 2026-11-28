<!-- https://docs.dseams.info/reference/nix.md -->

# Nix flake {#nix-flake}

The flake builds `libyodaLib` and the `seams` CLI with Meson. It does
not build Python or Lua. There are no overlays.

```bash
nix build                  # ./result/bin/seams
nix run . -- read input/traj/exampleTraj.lammpstrj
nix develop                # inputsFrom seams, plus gdb, ccache, clang-tools
nix flake update           # refresh the nixpkgs pin
nix fmt
```

`nix build` sets `doCheck`, so the build runs the Meson suite.
That suite includes the CLI tests `seams_read`, `seams_help`,
`seams_version`, `seams_features`, and `seams_print_config`.

## Flake outputs {#flake-outputs}

From `flake.nix`. Supported systems: `x86_64-linux`,
`aarch64-linux`, `x86_64-darwin`, `aarch64-darwin`.

| output | value |
| --- | --- |
| packages.<system>.seams | nix/package.nix derivation |
| packages.<system>.default | same as seams |
| apps.<system>.default | ${packages.default}/bin/seams |
| checks.<system>.seams | the seams package |
| checks.<system>.default | same as seams |
| devShells.<system>.default | inputsFrom the package; extra packages gdb, ccache, clang-tools |
| formatter.<system> | nixfmt |

The only flake input is `nixpkgs` (`github:NixOS/nixpkgs/nixos-unstable`).

## What the package enables {#what-the-package-enables}

From `nix/package.nix` (`pname = "seams"`, version set to the project version):

| on | Eigen, BLAS, LAPACK, Catch2 3, Highway, spglib, nauty, and the OpenBabel CSM program |
| --- | --- |
| engines | LAMMPS, GROMACS, OpenMM, and PLUMED from nixpkgs |
| neighbour geometry | vesin, linkcell, and minimage from fixed source archives |
| OpenMP | dependency('openmp') probe; Clang also gets llvmPackages.openmp |
| Meson flags | tests, CLI, packaged analysis and engine backends enabled; unavailable backends disabled |
| wrap mode | nodownload (Meson setup hook) |

The package uses the clean flake source. It stages checksum-pinned
trajectory data and binding documentation sources for the sandboxed checks.

The backward-cpp wrap uses a pinned source archive with its Meson overlay.

## What stays off {#what-stays-off}

| with_mpi | false |
| --- | --- |
| with_openmp_offload | false |
| readcon-core, readcon-db | no package on this nixpkgs pin |
| chemfiles, IRA, sphericart | no package on this nixpkgs pin |

`seams --features` on a flake build therefore reports Highway (and
OpenMP when it linked), spglib, nauty, vesin, linkcell, and minimage as
enabled. Chemfiles, readcon-core, IRA/SOFI, sphericart, and MPI stay
disabled because this nixpkgs pin has no package for them or the option
is intentionally off.

## Front-end flakes {#front-end-flakes}

[yodaStruct](https://github.com/d-SEAMS/yodaStruct) and
[PydSEAMSlib](https://github.com/d-SEAMS/PydSEAMSlib) each take this
repository as a source input (`flake = false`) and pass that tree
into their package. The yodaStruct flake builds it as a Meson
subproject.

## Cachix {#cachix}

The [dseams](https://dseams.cachix.org/) cache is optional (not a
flake input):

```bash
nix-env -iA cachix -f https://cachix.org/api/v1/install
cachix use dseams
```
