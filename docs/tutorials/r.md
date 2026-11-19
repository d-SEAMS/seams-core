# R

Classify ice Ih from R. The recording uses `input/traj/ice1h.lammpstrj`: `seams read` prints `nop 128`, and published CHILL+ at a 3.5 angstrom cutoff prints `nop 128 hexagonal 128`.

```{eval-rst}
.. include:: ../videos/embeds/r.rst
```

The binding is the R package `dseams` in [d-SEAMS/dseams-r](https://github.com/d-SEAMS/dseams-r). Install it from GitHub:

```text
remotes::install_github("d-SEAMS/dseams-r")
```

The package links `libseams_lang` (`bindings/cxx/seams_lang.h`).

The same two engine commands, without the package:

```bash
seams read input/traj/ice1h.lammpstrj
seams chill-plus input/traj/ice1h.lammpstrj --cutoff 3.5 --type 1
```

The package name is the {doc}`R API </stable/r>`.
