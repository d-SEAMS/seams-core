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

`build_tarballs.jl` sets `v2_14_0_tag_sha` to the commit v2.14.0 tags,
`acc8bbc803f55bac39493c85300f87a8ce302829` (`git rev-parse
v2.14.0^{commit}`). The Yggdrasil pull request carries that pin. Register
`libseams_jll` before Dseams: the General registry resolves
`libseams_jll` when Registrator runs.

## The patches

`bundled/patches` records the portability diff that is already contained
in `v2_14_0_tag_sha` (non-x86 hosts, musl, libc++, optional BLAS). The
recipe does not apply those patches. They stay in the tree as the record
of that diff.
