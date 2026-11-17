# Primitive rings

Shot list for a silent picture track. Draft narration is below and in narration.txt. Final wording, edit, and caption timing are recorded later.

Pixi environment: seams-core for the binary.

Input data: `input/traj/exampleTraj.lammpstrj`, `input/traj/ice1h.lammpstrj`.

Prompt on the terminal shots is `$ `. Commands use the build through `$SEAMS`, `$SEAMS_BUILD`, `$SEAMS_BUILD_SRC`, `$SEAMS_LIBDIR`, `$SEAMS_LANGDIR`, and `$VIDEO_OUT`.

Tests: `docs_video_tapes`.

Picture: 1920x1080, 30 fps, crf 16, one vhs scene file per shot so the lead-in can be retimed. Each shot starts 0.30 s before the draft line. Captions are a `<track>`, in the empty lower band (`line:95%`), and are not burned in. A later mix uses `-c:a aac -aac_coder fast`.

## Scene 1. exampleTraj census

Lead-in: 0.30 s. Segment length: 14.00 s.

Tests: `docs_video_tapes`.

Commands:

```bash
"$SEAMS" rings input/traj/exampleTraj.lammpstrj --type 2 --cutoff 3.5 --max-ring 6
```

Expected output (substrings the run must print):

```text
nop 250 rings 511
3 37
4 273
5 34
6 167
```

Draft narration:

exampleTraj.lammpstrj is 250 oxygens. seams rings prints nop 250 rings 511, then the census 3 37, 4 273, 5 34, and 6 167.

## Scene 2. Ice Ih six-rings

Lead-in: 0.30 s. Segment length: 12.00 s.

Tests: `docs_video_tapes`.

Commands:

```bash
"$SEAMS" rings input/traj/ice1h.lammpstrj --type 1 --cutoff 3.5 --max-ring 6
```

Expected output (substrings the run must print):

```text
nop 128 rings 256
6 256
```

Draft narration:

The ice Ih cell prints nop 128 rings 256, and the histogram line is 6 256. Those are the six-rings.

## End card

Goswami, R., Goswami, A., & Singh, J. K. (2020). d-SEAMS: Deferred Structural Elucidation Analysis for Molecular Simulations. Journal of Chemical Information and Modeling. https://doi.org/10.1021/acs.jcim.0c00031

d-SEAMS 2.0: Incremental ring and cage algorithms for ice structure analysis. Rohit Goswami, Amrita Goswami, Ruhila Goswami, Jayant K. Singh, to be submitted soon.

Questions, issues, discussions and pull requests welcome
