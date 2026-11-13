# The seams command on ice Ih

Shot list for a silent picture track. Draft narration is below and in narration.txt. Final wording, edit, and caption timing are recorded later.

Pixi environment: seams-core (`pixi run build`).

Input data: `input/traj/ice1h.lammpstrj`.

Prompt on the terminal shots is `$ `. Commands use the build through `$SEAMS`, `$SEAMS_BUILD`, `$SEAMS_BUILD_SRC`, `$SEAMS_LIBDIR`, `$SEAMS_LANGDIR`, and `$VIDEO_OUT`.

Tests: `docs_engine_openmm`, `docs_video_tapes`.

Picture: 1920x1080, 30 fps, crf 16, one vhs scene file per shot so the lead-in can be retimed. Each shot starts 0.20 to 0.45 s before the draft line. Captions are a `<track>`, in the empty lower band (`line:95%`), and are not burned in. A later mix uses `-c:a aac -aac_coder fast`.

## Scene 1. Install

Lead-in: 0.30 s. Segment length: 15.25 s.

Tests: `docs_video_tapes`.

Commands:

```bash
printf '%s\n' 'pixi run build' 'meson setup build' 'meson compile -C build'
```

Expected output (substrings the run must print):

```text
pixi run build
```

Draft narration:

The command-line tool is the seams binary from the seams-core pixi environment. pixi run build is meson setup build followed by meson compile -C build. After that, seams is the program on screen.

## Scene 2. Load ice Ih

Lead-in: 0.25 s. Segment length: 10.20 s.

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

I point seams read at input/traj/ice1h.lammpstrj. The line is nop 128 frame 1. One frame, 128 oxygens. docs_engine_openmm asserts this cell.

## Scene 3. Neighbours

Lead-in: 0.40 s. Segment length: 13.27 s.

Tests: `docs_video_tapes`.

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

seams cn with cutoff 3.5 and type 1 prints nI 128. Every oxygen has neighbours inside the cutoff. That list is the input to rings and to CHILL+.

## Scene 4. Rings

Lead-in: 0.20 s. Segment length: 9.73 s.

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

seams rings prints rings 256 and the histogram 6 256. The cell is six-rings only. docs_engine_openmm calls this six-rings 256.

## Scene 5. q6 and CHILL+

Lead-in: 0.35 s. Segment length: 11.13 s.

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

seams steinhardt for l 6 prints q6 0.625855. seams chill-plus prints hexagonal 128. The frame is entirely hexagonal ice under the published rule.

## Scene 6. Read the result

Lead-in: 0.45 s. Segment length: 17.07 s.

Tests: `docs_video_tapes`, `docs_engine_openmm`.

Commands:

```bash
seams chill-plus input/traj/ice1h.lammpstrj --cutoff 3.5 --type 1
seams rings input/traj/ice1h.lammpstrj --cutoff 3.5 --type 1
```

Expected output (substrings the run must print):

```text
hexagonal 128
6 256
```

Draft narration:

I read the four lines as one result. nop 128, nI 128, q6 0.625855, hexagonal 128, and 6 256. If a later change moves any of those, this ice Ih cell is the first place I look.
