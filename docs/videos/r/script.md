# R on ice Ih

Shot list for a silent picture track. Draft narration is below and in narration.txt. Final wording, edit, and caption timing are recorded later.

Pixi environment: seams-core (`pixi run build` in pixi.toml). The compiler and the language tool are on PATH.

Input data: `input/traj/ice1h.lammpstrj`, `input/traj/genice_sI.lammpstrj`.

Prompt on the terminal shots is `$ `. Commands use the build through `$SEAMS`, `$SEAMS_BUILD`, `$SEAMS_BUILD_SRC`, `$SEAMS_LIBDIR`, `$SEAMS_LANGDIR`, and `$VIDEO_OUT`.

Tests: `docs_engine_openmm`, `docs_video_tapes`.

Picture: 1920x1080, 30 fps, crf 16, one vhs scene file per shot so the lead-in can be retimed. Each shot starts 0.20 to 0.45 s before the draft line. Captions are a `<track>`, in the empty lower band (`line:95%`), and are not burned in. A later mix uses `-c:a aac -aac_coder fast`.

## Scene 1. Install

Lead-in: 0.30 s. Segment length: 24.00 s.

Tests: `docs_video_tapes`.

Commands:

```bash
printf '%s\n' 'pixi run build' 'meson setup build' 'meson compile -C build'
```

Expected output (substrings the run must print):

```text
pixi run build
meson setup build
meson compile -C build
```

Draft narration:

I build this tree with the seams-core pixi environment. The two lines on screen are the install: meson setup build, then meson compile -C build. That produces the seams command and libseams, which is the library the R call sits on. I keep that build and use it for every count in this recording.

## Scene 2. Load ice Ih

Lead-in: 0.25 s. Segment length: 15.62 s.

Tests: `docs_video_tapes`, `docs_engine_openmm`.

Commands:

```bash
seams read input/traj/ice1h.lammpstrj
```

Expected output (substrings the run must print):

```text
nop 128
frame 1
```

Draft narration:

The input is input/traj/ice1h.lammpstrj. seams read prints nop 128 frame 1. Those are the 128 oxygens in one frame. This is the cell docs_engine_openmm scores, and it is the cell I hand to R.

## Scene 3. Neighbours

Lead-in: 0.40 s. Segment length: 17.85 s.

Tests: `docs_video_tapes`, `docs_engine_openmm`.

Commands:

```bash
seams cn input/traj/ice1h.lammpstrj --cutoff 3.5 --type 1
```

Expected output (substrings the run must print):

```text
cutoff 3.5
nI 128
```

Draft narration:

Neighbours are a cutoff of 3.5 angstroms on type 1. seams cn prints nI 128, so every oxygen is in the list. Ice Ih keeps a four-neighbour shell inside that cutoff. The R neighbour call walks this same list.

## Scene 4. Rings

Lead-in: 0.20 s. Segment length: 15.57 s.

Tests: `docs_video_tapes`, `docs_engine_openmm`.

Commands:

```bash
seams rings input/traj/ice1h.lammpstrj --cutoff 3.5 --type 1
```

Expected output (substrings the run must print):

```text
rings 256
6 256
```

Draft narration:

Primitive rings use that list. The histogram is 6 256, and the summary is rings 256. Two six-rings for each of the 128 waters, and no other size. docs_engine_openmm names that count six-rings 256.

## Scene 5. q6 and CHILL+

Lead-in: 0.35 s. Segment length: 16.13 s.

Tests: `docs_video_tapes`, `docs_engine_openmm`.

Commands:

```bash
seams steinhardt input/traj/ice1h.lammpstrj --cutoff 3.5 --l 6 --type 1
seams chill-plus input/traj/ice1h.lammpstrj --cutoff 3.5 --type 1
```

Expected output (substrings the run must print):

```text
q6 0.625855
hexagonal 128
```

Draft narration:

Local q6 on the same cutoff prints q6 0.625855. Published CHILL+ then prints hexagonal 128. Every water is hexagonal. I want the R result to show the same two figures before I trust the wrapper.

## Scene 6. Read the result

Lead-in: 0.45 s. Segment length: 16.65 s.

Tests: `docs_video_tapes`, `docs_engine_openmm`.

Commands:

```bash
Rscript docs/examples/minimal_r.R input/traj/genice_sI.lammpstrj
```

Expected output (substrings the run must print):

```text
r example ok
```

Draft narration:

The R example is docs/examples/minimal_r.R on input/traj/genice_sI.lammpstrj, with the dseams package installed from bindings/r/dseams. It prints r example ok. The ice Ih census from seams is still nop 128, q6 0.625855, hexagonal 128, and rings 256.

## End card

Goswami, R., Goswami, A., & Singh, J. K. (2020). d-SEAMS: Deferred Structural Elucidation Analysis for Molecular Simulations. Journal of Chemical Information and Modeling. https://doi.org/10.1021/acs.jcim.0c00031

d-SEAMS 2.0: Incremental ring and cage algorithms for ice structure analysis. Rohit Goswami, Amrita Goswami, Ruhila Goswami, Jayant K. Singh, to be submitted soon.

Questions, issues, discussions and pull requests welcome
