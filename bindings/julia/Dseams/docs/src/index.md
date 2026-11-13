# Dseams

`Dseams` calls `libseams_lang`. Set `SEAMS_LANG_LIB` to that shared library.
`SEAMS_LIB` is `libseams` for the `seams_*` calls included from the C API.

```julia
using Dseams
frame = read_frame("ice.lammpstrj"; type = 1)
ice = chill_plus(frame)
println(ice)
println(ice.counts["cubic"])
```

A `Frame` finalizes itself. A failing call throws `SeamsError`.
The message is the string stored for that call.

## Reference

```@autodocs
Modules = [Dseams]
```
