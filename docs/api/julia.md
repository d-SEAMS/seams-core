# Julia

| package | calls | library |
| --- | --- | --- |
| `Dseams` | the language ABI | `libseams_lang` (`bindings/cxx/seams_lang.h`) |
| `Seams` | the C API | `libseams` (`include/seams.h`) |

Both live in [d-SEAMS/Dseams.jl](https://github.com/d-SEAMS/Dseams.jl). Neither is in the General registry. The README's local-engine section names `SEAMS_LANG_LIB` and `SEAMS_LIB`.

The first classification is the {doc}`Julia tutorial </fronts/julia>`.
