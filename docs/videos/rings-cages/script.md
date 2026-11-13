# Rings and sI, sII cages

Shot list for a silent picture track. Draft narration is below and in narration.txt. Final wording, edit, and caption timing are recorded later.

Pixi environment: seams-core for the binary, and the dseams-engines manifest (openmm, genice2) for the tutorial driver.

Input data: `docs/tutorials/paper_science.py rings-full and cages-full`.

Prompt on the terminal shots is `$ `. Commands use the build through `$SEAMS`, `$SEAMS_BUILD`, `$SEAMS_BUILD_SRC`, `$SEAMS_LIBDIR`, `$SEAMS_LANGDIR`, and `$VIDEO_OUT`.

Tests: `docs_paper_cages_full`, `docs_paper_rings_full`, `docs_video_tapes`.

Picture: 1920x1080, 30 fps, crf 16, one vhs scene file per shot so the lead-in can be retimed. Each shot starts 0.30 to 0.35 s before the draft line. Captions are a `<track>`, in the empty lower band (`line:95%`), and are not burned in. A later mix uses `-c:a aac -aac_coder fast`.

## Scene 1. Ice Ih rings

Lead-in: 0.30 s. Segment length: 12.33 s.

Tests: `docs_video_tapes`, `docs_paper_rings_full`.

Commands:

```bash
python3 docs/tutorials/paper_science.py rings-full "$SEAMS" "$SEAMS_LIBDIR" "$VIDEO_OUT/rings"
```

Expected output (substrings the run must print):

```text
rings full nop 4096 rings 8192
rings full histogram 6 8192
```

Draft narration:

The full ice Ih cell is 4096 molecules. rings-full prints rings full nop 4096 rings 8192 and rings full histogram 6 8192. docs_paper_rings_full asserts both lines.

## Scene 2. sI and sII cages

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
