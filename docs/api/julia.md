# Julia

| package | calls | library |
| --- | --- | --- |
| `Dseams` | the language ABI | `libseams_lang` (`bindings/cxx/seams_lang.h`) |
| `Seams` | the C API | `libseams` (`include/seams.h`) |

Both live in [d-SEAMS/Dseams.jl](https://github.com/d-SEAMS/Dseams.jl). Install as a Git dependency:

```bash
julia -e 'using Pkg; Pkg.add(url="https://github.com/d-SEAMS/Dseams.jl")'
```

The first classification is the {doc}`Julia tutorial </fronts/julia>`.
