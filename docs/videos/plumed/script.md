# PLUMED driver

Shot list for a silent picture track. Draft narration is below and in narration.txt. Final wording, edit, and caption timing are recorded later.

Pixi environment: dseams-engines (lammps, gromacs, plumed 2.10) plus seams-core for the hook.

Input data: `input/traj/ice1h.lammpstrj`, `engines/plumed/plumed.dat`.

Prompt on the terminal shots is `$ `. Commands use the build through `$SEAMS`, `$SEAMS_BUILD`, `$SEAMS_BUILD_SRC`, `$SEAMS_LIBDIR`, `$SEAMS_LANGDIR`, and `$VIDEO_OUT`.

Tests: `docs_engine_plumed`, `docs_video_tapes`, `plumed_driver`.

Picture: 1920x1080, 30 fps, crf 16, one vhs scene file per shot so the lead-in can be retimed. Each shot starts 0.32 to 0.32 s before the draft line. Captions are a `<track>`, in the empty lower band (`line:95%`), and are not burned in. A later mix uses `-c:a aac -aac_coder fast`.

## Scene 1. Driver on ice Ih

Lead-in: 0.32 s. Segment length: 16.10 s.

Tests: `docs_video_tapes`, `docs_engine_plumed`, `plumed_driver`.

Commands:

```bash
PLUMED_NUM_THREADS=1 python3 docs/howto/test_engines.py plumed . "$SEAMS_BUILD"
```

Expected output (substrings the run must print):

```text
plumed hexagonal 128 six-rings 256 q6 0.625855
```

Draft narration:

PLUMED_NUM_THREADS is set to one thread before plumed driver starts. docs_engine_plumed loads libSeamsPlumed.so, writes an xyz in nanometres, and prints plumed hexagonal 128 six-rings 256 q6 0.625855. plumed_driver checks the cubic cell the same way.
