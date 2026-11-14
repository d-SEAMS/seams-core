# Primitive rings

Shot list for a silent picture track. Draft narration is below and in narration.txt. Final wording, edit, and caption timing are recorded later.

Pixi environment: seams-core for the binary, and the dseams-engines manifest (openmm, genice2) for the tutorial driver.

Input data: `docs/tutorials/paper_science.py rings-full`.

Prompt on the terminal shots is `$ `. Commands use the build through `$SEAMS`, `$SEAMS_BUILD`, `$SEAMS_BUILD_SRC`, `$SEAMS_LIBDIR`, `$SEAMS_LANGDIR`, and `$VIDEO_OUT`.

Tests: `docs_paper_rings_full`, `docs_video_tapes`.

Picture: 1920x1080, 30 fps, crf 16, one vhs scene file per shot so the lead-in can be retimed. Each shot starts 0.30 s before the draft line. Captions are a `<track>`, in the empty lower band (`line:95%`), and are not burned in. A later mix uses `-c:a aac -aac_coder fast`.

## Scene 1. Ice Ih rings

Lead-in: 0.30 s. Segment length: 10.67 s.

Tests: `docs_video_tapes`, `docs_paper_rings_full`.

Commands:

```bash
python3 docs/tutorials/paper_science.py rings-full "$SEAMS" "$SEAMS_LIBDIR" "$VIDEO_OUT/rings"
```

Expected output (substrings the run must print):

```text
rings full input molecules 4096 frames 1 file-frames 11 file input/traj/mW_cubic.lammpstrj
rings full paper molecules 4096 frames 1 rings 8192
rings full nop 4096 rings 8192
rings full histogram 6 8192
```

Draft narration:

The paper ring frame is input/traj/mW_cubic.lammpstrj. It holds 11 frames of 4096 molecules. Frame 1 has 8192 six-rings. docs_paper_rings_full asserts those lines.

## End card

Goswami, R., Goswami, A., & Singh, J. K. (2020). d-SEAMS: Deferred Structural Elucidation Analysis for Molecular Simulations. Journal of Chemical Information and Modeling. https://doi.org/10.1021/acs.jcim.0c00031

d-SEAMS 2.0: Incremental ring and cage algorithms for ice structure analysis. Rohit Goswami, Amrita Goswami, Ruhila Goswami, Jayant K. Singh, to be submitted soon.
