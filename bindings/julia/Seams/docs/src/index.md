# Seams

`Seams` calls `libseams` (`include/seams.h`). Set `SEAMS_LIB` to
`libseams.so`.

```julia
using Seams
set_num_threads(1)
sys = System()
set_ortho!(sys, xyz, (lx, ly, lz))
labels, counts = chill_plus(sys; rule = PUBLISHED)
q4, q6 = q4q6(sys)
six, fraction = six_rings(sys)
```

`System` is finalized with `seams_system_destroy`.
`ChillLabels`, `Q4Q6Pair`, and `SixRings` iterate in field order.
A failing call throws `SeamsError`, whose message is `seams_handle_error`.

## Reference

```@autodocs
Modules = [Seams]
```
