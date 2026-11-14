# Migrate from 1.x

Shot list for a silent picture track. Draft narration is below and in narration.txt. Final wording, edit, and caption timing are recorded later.

Pixi environment: seams-core (`pixi run build`).

Input data: `input/traj/exampleTraj.lammpstrj`.

Prompt on the terminal shots is `$ `. Commands use the build through `$SEAMS`, `$SEAMS_BUILD`, `$SEAMS_BUILD_SRC`, `$SEAMS_LIBDIR`, `$SEAMS_LANGDIR`, and `$VIDEO_OUT`.

Tests: `docs_migrate`, `docs_video_tapes`.

Picture: 1920x1080, 30 fps, crf 16, one vhs scene file per shot so the lead-in can be retimed. Each shot starts 0.30 to 0.35 s before the draft line. Captions are a `<track>`, in the empty lower band (`line:95%`), and are not burned in. A later mix uses `-c:a aac -aac_coder fast`.

## Scene 1. Hydrogen bonds and coordination

Lead-in: 0.30 s. Segment length: 13.17 s.

Tests: `docs_video_tapes`, `docs_migrate`.

Commands:

```bash
seams hbonds input/traj/exampleTraj.lammpstrj --type 2 --htype 1
seams cn input/traj/exampleTraj.lammpstrj --pair-types 1,2 --cutoff 4.5
```

Expected output (substrings the run must print):

```text
nop 250 hbonds 460
cn 9.144
cutoff 4.5
```

Draft narration:

The migration check uses input/traj/exampleTraj.lammpstrj. seams hbonds prints nop 250 hbonds 460. seams cn for the 1,2 pair at cutoff 4.5 prints cn 9.144. docs_migrate requires both lines.

## Scene 2. Fingerprint, pairs, domains

Lead-in: 0.35 s. Segment length: 20.30 s.

Tests: `docs_video_tapes`, `docs_migrate`.

Commands:

```bash
seams fingerprint input/traj/exampleTraj.lammpstrj --type 2
seams pairs input/traj/exampleTraj.lammpstrj --site 1=cationHead,2=anion
seams domains input/traj/exampleTraj.lammpstrj --site 1=cationHead,2=anion --subset polar
```

Expected output (substrings the run must print):

```text
classes 89
count 0
largest 750
```

Draft narration:

seams fingerprint prints classes 89. seams pairs prints count 0. seams domains on the polar subset prints largest 750. docs_migrate requires classes 89, count 0, and largest 750. The old name and the current command in the guide are the two sides of that check.

## End card

Goswami, R., Goswami, A., & Singh, J. K. (2020). d-SEAMS: Deferred Structural Elucidation Analysis for Molecular Simulations. Journal of Chemical Information and Modeling. https://doi.org/10.1021/acs.jcim.0c00031

d-SEAMS 2.0: Incremental ring and cage algorithms for ice structure analysis. Rohit Goswami, Amrita Goswami, Ruhila Goswami, Jayant K. Singh, to be submitted soon.
