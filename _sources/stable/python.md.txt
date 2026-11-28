# Python

`pydseamslib` provides `import dseams` and `import pydseams`.

| name | role | library |
| --- | --- | --- |
| `dseams` | array calls | `libseams` (`include/seams.h`) |
| `pydseams` | Frame API | `libseams` and `libseams_lang` (`bindings/cxx/seams_lang.h`) |

The package and its symbol list are [d-SEAMS/PydSEAMSlib](https://github.com/d-SEAMS/PydSEAMSlib). The rendered reference is [PydSEAMSlib](https://d-seams.github.io/PydSEAMSlib/).

```bash
pip install pydseamslib
```

`dseams.chill_plus(coordinates, box)` returns one label per atom and nine counts. Index 0 is cubic and index 1 is hexagonal. On the eight-molecule cubic cell in the quickstart, index 0 is 8. On `input/traj/ice1h.lammpstrj` (128 oxygens, LAMMPS type 1) the same call returns 128 at index 1, and the CLI prints:

```
nop 128 frame 1 box 15.6457 14.7071 18.073
nop 128 hexagonal 128
```

The first classification is the {doc}`Python tutorial </fronts/python>`.
