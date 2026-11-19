# Julia

| package | calls | library |
| --- | --- | --- |
| `Dseams` | the language ABI | `libseams_lang` (`bindings/cxx/seams_lang.h`) |
| `Seams` | the C API | `libseams` (`include/seams.h`) |

Both live in [d-SEAMS/Dseams.jl](https://github.com/d-SEAMS/Dseams.jl). Install as a Git dependency:

```bash
julia -e 'using Pkg; Pkg.add(url="https://github.com/d-SEAMS/Dseams.jl")'
```

`Dseams` reads a frame through `libseams_lang`. On `input/traj/ice1h.lammpstrj` (128 oxygens, type 1) published CHILL+ at 3.5 angstrom is 128 hexagonal:

```
nop 128 frame 1 box 15.6457 14.7071 18.073
nop 128 hexagonal 128
```

The first classification is the {doc}`Julia tutorial </fronts/julia>`.
