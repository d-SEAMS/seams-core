# Julia

The Julia binding is the Julia package `Dseams`. It lives in its own repository,
[d-SEAMS/Dseams.jl](https://github.com/d-SEAMS/Dseams.jl), with its API reference, tutorials,
tests and release notes. It links the engine through the C API
(`libseams`, `include/seams.h`) and the C ABI `libseams_lang`
(`bindings/cxx/seams_lang.h`) that this repository builds and installs.

Install the current source as a Git dependency:

```bash
julia -e 'using Pkg; Pkg.add(url="https://github.com/d-SEAMS/Dseams.jl")'
```
