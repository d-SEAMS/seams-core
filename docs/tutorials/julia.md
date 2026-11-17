# Julia

Classify ice Ih from Julia. The recording uses `input/traj/ice1h.lammpstrj`: `seams read` prints `nop 128`, and published CHILL+ at a 3.5 angstrom cutoff prints `nop 128 hexagonal 128`.

```{eval-rst}
.. include:: ../videos/embeds/julia.rst
```

The binding is the Julia package `Dseams` in [d-SEAMS/Dseams.jl](https://github.com/d-SEAMS/Dseams.jl). It is not in the General registry. The README's local-engine section is the install, and its first example classifies the bundled ice Ih frame as 128 hexagonal.

`Dseams` calls `libseams_lang` (`bindings/cxx/seams_lang.h`). The same repository's `Seams` package calls `libseams` (`include/seams.h`).

The same two engine commands, without the binding:

```bash
seams read input/traj/ice1h.lammpstrj
seams chill-plus input/traj/ice1h.lammpstrj --cutoff 3.5 --type 1
```

Package names are the {doc}`Julia API </stable/julia>`.
