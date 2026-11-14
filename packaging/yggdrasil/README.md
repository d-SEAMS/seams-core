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

## The patches

`bundled/patches` holds the seams-core changes made after the pinned
commit that the cross builds need (non-x86 hosts, musl, libc++, optional
BLAS). Once a seams-core release contains them, point the `GitSource` at
that release's commit and delete the patches.
