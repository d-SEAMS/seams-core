# Python

Classify ice Ih from Python. The recording uses `input/traj/ice1h.lammpstrj`: `seams read` prints `nop 128`, and published CHILL+ at a 3.5 angstrom cutoff prints `nop 128 hexagonal 128`.

```{eval-rst}
.. include:: ../videos/embeds/python.rst
```

The binding is [d-SEAMS/PydSEAMSlib](https://github.com/d-SEAMS/PydSEAMSlib). `import dseams` is the array API over `libseams` (`include/seams.h`). `import pydseams` is the Frame API. The worked example and the symbol list are the [PydSEAMSlib documentation](https://d-seams.github.io/PydSEAMSlib/).

```bash
pip install 'pydseamslib==2.14.0'
```

The same two engine commands, without the binding, print the ice Ih census. The box is orthorhombic and every oxygen is hexagonal:

```bash
seams read input/traj/ice1h.lammpstrj
seams chill-plus input/traj/ice1h.lammpstrj --cutoff 3.5 --type 1
```

```
nop 128 frame 1 box 15.6457 14.7071 18.073
nop 128 hexagonal 128
```

Package names and the libraries they link are the {doc}`Python API </stable/python>`.
