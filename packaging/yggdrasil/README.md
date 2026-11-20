# libseams_jll

`L/libseams` is the Yggdrasil recipe for `libseams_jll`, the engine
libraries (`libyodaLib`, `libseams`, `libseams_lang`) and the `seams`
CLI that the Julia package Dseams loads. It mirrors the Yggdrasil layout:
submitting is copying `L/libseams` into a Yggdrasil checkout.

## Building it here

The recipe includes Yggdrasil's `platforms/macos_sdks.jl`, two levels up
from the recipe. Outside Yggdrasil, fetch it next to `L/` first:

```sh
mkdir -p packaging/yggdrasil/platforms
curl -fsSL -o packaging/yggdrasil/platforms/macos_sdks.jl \
  https://raw.githubusercontent.com/JuliaPackaging/Yggdrasil/master/platforms/macos_sdks.jl
cd packaging/yggdrasil/L/libseams
julia --project=@bb -e 'using Pkg; Pkg.add("BinaryBuilder")'
julia --project=@bb build_tarballs.jl --verbose --deploy=local x86_64-linux-gnu-cxx11
```

The audit step of the BinaryBuilder 0.6.6 release has no riscv64 entry and
stops on that platform after it builds; Yggdrasil's CI runs a newer
BinaryBuilder. Locally, list the other platforms explicitly.

`--deploy=local` writes the JLL to `~/.julia/dev/libseams_jll`. To test
Dseams against it:

```sh
julia --project=ENV bindings/julia/Dseams/dev/develop_jll.jl ~/.julia/dev/libseams_jll
julia --project=ENV -e 'using Pkg; Pkg.develop(path="bindings/julia/Dseams"); Pkg.test("Dseams")'
```

## The source pin

`build_tarballs.jl` sets `tag_sha` as the stand-in for the `v<version>`
tag matching `version` in that file (currently the release-candidate
commit `c1efa8c734214e9117cb77e89258bc9fdd214a6c`). After creating the
release tag (`VERSION=$(python3 scripts/check_release_version.py --print)`),
replace the constant with `git rev-parse "v${VERSION}^{commit}"` when the
tag points somewhere else, then open the Yggdrasil pull request. Register
`libseams_jll` before Dseams: the General registry resolves
`libseams_jll` when Registrator runs.

## The patches

`bundled/patches` records the portability diff that is already contained
in `tag_sha` (non-x86 hosts, musl, libc++, optional BLAS). The
recipe does not apply those patches. They stay in the tree as the record
of that diff.
