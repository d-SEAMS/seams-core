# GROMACS gmx-seams

Shot list for a silent picture track. Draft narration is below and in narration.txt. Final wording, edit, and caption timing are recorded later.

Pixi environment: dseams-engines (lammps, gromacs, plumed 2.10) plus seams-core for the hook.

Input data: `engines/gromacs/ice.gro`, classified by gmx-seams after gmx mdrun.

Prompt on the terminal shots is `$ `. Commands use the build through `$SEAMS`, `$SEAMS_BUILD`, `$SEAMS_BUILD_SRC`, `$SEAMS_LIBDIR`, `$SEAMS_LANGDIR`, and `$VIDEO_OUT`.

Tests: `docs_engine_gromacs`, `docs_video_tapes`, `gromacs_mdrun`.

Picture: 1920x1080, 30 fps, crf 16, one vhs scene file per shot so the lead-in can be retimed. Each shot starts 0.30 to 0.30 s before the draft line. Captions are a `<track>`, in the empty lower band (`line:95%`), and are not burned in. A later mix uses `-c:a aac -aac_coder fast`.

## Scene 1. mdrun, then gmx-seams

Lead-in: 0.30 s. Segment length: 22.00 s.

Tests: `docs_video_tapes`, `gromacs_mdrun`.

Commands:

```bash
engines/gromacs/run_mdrun.sh "$SEAMS_BUILD/engines/gmx-seams" "$SEAMS_BUILD_SRC"
```

Expected output (substrings the run must print):

```text
gromacs frames 5 ice 320 atoms 320
gromacs hexagonal 0 six-rings 128 q6 0.628519
```

Draft narration:

engines/gromacs/run_mdrun.sh runs gmx grompp, then gmx mdrun for 4 steps, then gmx-seams on the xtc mdrun wrote. The census line is gromacs frames 5 ice 320 atoms 320. ice 320 is the cubic total over 5 frames of 64. The last frame is gromacs hexagonal 0 six-rings 128 q6 0.628519. hexagonal, six-rings, and q6 are last-frame values. gmx-seams is post-processing after mdrun.

## End card

Goswami, R., Goswami, A., & Singh, J. K. (2020). d-SEAMS: Deferred Structural Elucidation Analysis for Molecular Simulations. Journal of Chemical Information and Modeling. https://doi.org/10.1021/acs.jcim.0c00031

d-SEAMS 2.0: Incremental ring and cage algorithms for ice structure analysis. Rohit Goswami, Amrita Goswami, Ruhila Goswami, Jayant K. Singh, to be submitted soon.

Questions, issues, discussions and pull requests welcome
