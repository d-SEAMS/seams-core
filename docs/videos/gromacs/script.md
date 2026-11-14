# GROMACS gmx-seams

Shot list for a silent picture track. Draft narration is below and in narration.txt. Final wording, edit, and caption timing are recorded later.

Pixi environment: dseams-engines (lammps, gromacs, plumed 2.10) plus seams-core for the hook.

Input data: `input/traj/ice1h.lammpstrj`, `engines/gromacs/ice.gro`.

Prompt on the terminal shots is `$ `. Commands use the build through `$SEAMS`, `$SEAMS_BUILD`, `$SEAMS_BUILD_SRC`, `$SEAMS_LIBDIR`, `$SEAMS_LANGDIR`, and `$VIDEO_OUT`.

Tests: `docs_engine_gromacs`, `docs_video_tapes`, `gromacs_mdrun`.

Picture: 1920x1080, 30 fps, crf 16, one vhs scene file per shot so the lead-in can be retimed. Each shot starts 0.30 to 0.30 s before the draft line. Captions are a `<track>`, in the empty lower band (`line:95%`), and are not burned in. A later mix uses `-c:a aac -aac_coder fast`.

## Scene 1. Ice Ih gro

Lead-in: 0.30 s. Segment length: 17.33 s.

Tests: `docs_video_tapes`, `docs_engine_gromacs`, `gromacs_mdrun`.

Commands:

```bash
python3 docs/howto/test_engines.py gromacs . "$SEAMS_BUILD"
```

Expected output (substrings the run must print):

```text
gromacs hexagonal 128 six-rings 256 q6 0.625855
```

Draft narration:

docs_engine_gromacs writes a one-frame gro from the ice Ih dump and runs gmx-seams with cutoff 3.5. The census line is gromacs hexagonal 128 six-rings 256 q6 0.625855. gromacs_mdrun then takes a short frozen trajectory through the same tool.

## End card

Goswami, R., Goswami, A., & Singh, J. K. (2020). d-SEAMS: Deferred Structural Elucidation Analysis for Molecular Simulations. Journal of Chemical Information and Modeling. https://doi.org/10.1021/acs.jcim.0c00031

d-SEAMS 2.0: Incremental ring and cage algorithms for ice structure analysis. Rohit Goswami, Amrita Goswami, Ruhila Goswami, Jayant K. Singh, to be submitted soon.
