# Fortran

Classify ice Ih from Fortran. The recording uses `input/traj/ice1h.lammpstrj`: `seams read` prints `nop 128`, and published CHILL+ at a 3.5 angstrom cutoff prints `nop 128 hexagonal 128`.

```{eval-rst}
.. include:: ../videos/embeds/fortran.rst
```

The binding is the fpm Git dependency `dseams` in [d-SEAMS/dseams-fortran](https://github.com/d-SEAMS/dseams-fortran). Programs use the module `seams`. `seams_c` is the `iso_c_binding` image of `include/seams.h`. The package links an installed `libseams`.

```toml
[dependencies]
dseams = { git = "https://github.com/d-SEAMS/dseams-fortran" }
```

The same two engine commands, without the package:

```bash
seams read input/traj/ice1h.lammpstrj
seams chill-plus input/traj/ice1h.lammpstrj --cutoff 3.5 --type 1
```

Module names are the {doc}`Fortran API </stable/fortran>`.
