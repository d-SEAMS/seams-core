# LAMMPS run_hook

Shot list for a silent picture track. Draft narration is below and in narration.txt. Final wording, edit, and caption timing are recorded later.

Pixi environment: dseams-engines (lammps, gromacs, plumed 2.10) plus seams-core for the hook.

Input data: `engines/lammps/in.dseams`, `engines/lammps/ice.data`, `input/traj/ice1h.lammpstrj`.

Prompt on the terminal shots is `$ `. Commands use the build through `$SEAMS`, `$SEAMS_BUILD`, `$SEAMS_BUILD_SRC`, `$SEAMS_LIBDIR`, `$SEAMS_LANGDIR`, and `$VIDEO_OUT`.

Tests: `docs_engine_lammps`, `docs_tutorial_insitu`, `docs_video_tapes`.

Picture: 1920x1080, 30 fps, crf 16, one vhs scene file per shot so the lead-in can be retimed. Each shot starts 0.30 to 0.35 s before the draft line. Captions are a `<track>`, in the empty lower band (`line:95%`), and are not burned in. A later mix uses `-c:a aac -aac_coder fast`.

## Scene 1. Cubic cell

Lead-in: 0.30 s. Segment length: 19.83 s.

Tests: `docs_video_tapes`, `docs_tutorial_insitu`.

Commands:

```bash
cd engines/lammps && "$SEAMS_BUILD/engines/run_hook" in.dseams
```

Expected output (substrings the run must print):

```text
lammps calls 5 stride 2 ice 320
lammps hexagonal 0 six-rings 128 q6 0.628539
```

Draft narration:

From engines/lammps the driver is the build's engines/run_hook and the deck is in.dseams. The cell has 64 atoms. The hook prints lammps calls 5 stride 2 ice 320, then lammps hexagonal 0 six-rings 128 q6 0.628539. docs_tutorial_insitu requires ice equal to calls times 64.

## Scene 2. Ice Ih census

Lead-in: 0.35 s. Segment length: 13.63 s.

Tests: `docs_video_tapes`, `docs_engine_lammps`.

Commands:

```bash
python3 docs/howto/test_engines.py lammps . "$SEAMS_BUILD"
```

Expected output (substrings the run must print):

```text
lammps hexagonal 128 six-rings 256 q6 0.625855
```

Draft narration:

docs_engine_lammps writes a data file from the ice Ih dump and runs the same hook for one step. The census line is lammps hexagonal 128 six-rings 256 q6 0.625855.
