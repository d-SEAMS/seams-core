# LAMMPS, GROMACS, OpenMM, and PLUMED

Shot list for a silent picture track. Draft narration is below and in narration.txt. Final wording, edit, and caption timing are recorded later.

Pixi environment: dseams-engines (openmm, and gromacs, lammps, plumed when those libraries are built). seams-core supplies libseams.

Input data: `input/traj/ice1h.lammpstrj`.

Prompt on the terminal shots is `$ `. Commands use the build through `$SEAMS`, `$SEAMS_BUILD`, `$SEAMS_BUILD_SRC`, `$SEAMS_LIBDIR`, `$SEAMS_LANGDIR`, and `$VIDEO_OUT`.

Tests: `docs_engine_openmm`, `docs_video_tapes`.

Picture: 1920x1080, 30 fps, crf 16, one vhs scene file per shot so the lead-in can be retimed. Each shot starts 0.30 to 0.40 s before the draft line. Captions are a `<track>`, in the empty lower band (`line:95%`), and are not burned in. A later mix uses `-c:a aac -aac_coder fast`.

## Scene 1. Shared ice Ih census

Lead-in: 0.30 s. Segment length: 19.83 s.

Tests: `docs_video_tapes`, `docs_engine_openmm`.

Commands:

```bash
python3 docs/howto/test_engines.py openmm . "$SEAMS_BUILD"
```

Expected output (substrings the run must print):

```text
openmm hexagonal 128 six-rings 256 q6 0.625855
```

Draft narration:

The four engine how-tos score the ice Ih dump. The locked line is hexagonal 128, six-rings 256, q6 0.625855. docs_engine_openmm asserts that line. On this build the OpenMM how-to is the one that runs, and it prints openmm hexagonal 128 six-rings 256 q6 0.625855.

## Scene 2. In-situ frame

Lead-in: 0.40 s. Segment length: 14.93 s.

Tests: `docs_video_tapes`, `docs_engine_openmm`.

Commands:

```bash
python3 docs/howto/test_engines.py insitu . "$SEAMS_BUILD"
```

Expected output (substrings the run must print):

```text
insitu hexagonal 128 six-rings 256
```

Draft narration:

The in-situ call reads the same dump through libseams_lang and prints insitu hexagonal 128 six-rings 256. LAMMPS, GROMACS, and PLUMED each have a shot list that runs the same ice Ih census.

## End card

Goswami, R., Goswami, A., & Singh, J. K. (2020). d-SEAMS: Deferred Structural Elucidation Analysis for Molecular Simulations. Journal of Chemical Information and Modeling. https://doi.org/10.1021/acs.jcim.0c00031

d-SEAMS 2.0: Incremental ring and cage algorithms for ice structure analysis. Rohit Goswami, Amrita Goswami, Ruhila Goswami, Jayant K. Singh, to be submitted soon.
