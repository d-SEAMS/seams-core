# sI, sII, and sH cages

Shot list for a silent picture track. Draft narration is below and in narration.txt. Final wording, edit, and caption timing are recorded later.

Pixi environment: seams-core for the binary, and the dseams-engines manifest (openmm, genice2) for the tutorial driver.

Input data: `docs/tutorials/paper_science.py cages-full`.

Prompt on the terminal shots is `$ `. Commands use the build through `$SEAMS`, `$SEAMS_BUILD`, `$SEAMS_BUILD_SRC`, `$SEAMS_LIBDIR`, `$SEAMS_LANGDIR`, and `$VIDEO_OUT`.

Tests: `docs_paper_cages_full`, `docs_video_tapes`.

Picture: 1920x1080, 30 fps, crf 16, one vhs scene file per shot so the lead-in can be retimed. Each shot starts 0.35 s before the draft line. Captions are a `<track>`, in the empty lower band (`line:95%`), and are not burned in. A later mix uses `-c:a aac -aac_coder fast`.

## Scene 1. Ice Ih hexagonal cages

Lead-in: 0.35 s. Segment length: 16.00 s.

Tests: `docs_video_tapes`.

Commands:

```bash
"$SEAMS" cages input/traj/ice1h.lammpstrj --graph seeded --signature hc+ddc
```

Expected output (substrings the run must print):

```text
nop 128 graph seeded signature 4:6,6:2 cages 64 atoms 128
nop 128 graph seeded signature 6:7 cages 0 atoms 0
```

Draft narration:

Ice Ih with the seeded hc+ddc signature prints nop 128 graph seeded signature 4:6,6:2 cages 64 atoms 128, then nop 128 graph seeded signature 6:7 cages 0 atoms 0. The first line is 64 hexagonal cages. The second search finds no double-diamond cages.

## Scene 2. Cubic double-diamond cages

Lead-in: 0.30 s. Segment length: 16.00 s.

Tests: `docs_video_tapes`.

Commands:

```bash
"$SEAMS" cages input/traj/mW_cubic.lammpstrj --type 1 --signature hc+ddc
```

Expected output (substrings the run must print):

```text
nop 4096 graph seeded signature 4:6,6:2 cages 0 atoms 0
nop 4096 graph seeded signature 6:7 cages 2048 atoms 4096
```

Draft narration:

The cubic cell prints nop 4096 graph seeded signature 4:6,6:2 cages 0 atoms 0, then nop 4096 graph seeded signature 6:7 cages 2048 atoms 4096. Those 2048 cages are the double-diamond cages.

## Scene 3. Clathrate cells

Lead-in: 0.35 s. Segment length: 15.30 s.

Tests: `docs_video_tapes`, `docs_paper_cages_full`.

Commands:

```bash
python3 docs/tutorials/paper_science.py cages-full "$SEAMS" "$SEAMS_LIBDIR" "$VIDEO_OUT/cages"
```

Expected output (substrings the run must print):

```text
cages full sI 16 48
cages full sII 64 32
cages full sH 24 16 8
```

Draft narration:

cages-full prints cages full sI 16 48, cages full sII 64 32, and cages full sH 24 16 8. The first two are the sI and sII clusters. docs_paper_cages_full asserts all three lines.

## End card

Goswami, R., Goswami, A., & Singh, J. K. (2020). d-SEAMS: Deferred Structural Elucidation Analysis for Molecular Simulations. Journal of Chemical Information and Modeling. https://doi.org/10.1021/acs.jcim.0c00031

d-SEAMS 2.0: Incremental ring and cage algorithms for ice structure analysis. Rohit Goswami, Amrita Goswami, Ruhila Goswami, Jayant K. Singh, to be submitted soon.

Questions, issues, discussions and pull requests welcome
