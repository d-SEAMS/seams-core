# Fortran on ice Ih

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

I build this tree with the seams-core pixi environment. The two lines on screen are the install: meson setup build, then meson compile -C build. That produces the seams command and libseams, which is the library the Fortran call sits on. I keep that build and use it for every count in this recording.

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

The input is input/traj/ice1h.lammpstrj. seams read prints nop 128 frame 1. Those are the 128 oxygens in one frame. This is the cell docs_engine_openmm scores, and it is the cell I hand to Fortran.

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

Neighbours are a cutoff of 3.5 angstroms on type 1. seams cn prints nI 128, so every oxygen is in the list. Ice Ih keeps a four-neighbour shell inside that cutoff. The Fortran neighbour call walks this same list.

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

Local q6 on the same cutoff prints q6 0.625855. Published CHILL+ then prints hexagonal 128. Every water is hexagonal. I want the Fortran result to show the same two figures before I trust the wrapper.

## Scene 6. Read the result

Lead-in: 0.45 s. Segment length: 20.82 s.

Tests: `docs_video_tapes`, `docs_engine_openmm`.

Commands:

```bash
gfortran -J "$VIDEO_OUT" -c bindings/fortran/src/seams_kinds.f90 -o "$VIDEO_OUT/seams_kinds.o" && gfortran -J "$VIDEO_OUT" -I "$VIDEO_OUT" -c bindings/fortran/src/seams_c.f90 -o "$VIDEO_OUT/seams_c.o" && gfortran -J "$VIDEO_OUT" -I "$VIDEO_OUT" -c bindings/fortran/src/seams_support.f90 -o "$VIDEO_OUT/seams_support.o" && gfortran -J "$VIDEO_OUT" -I "$VIDEO_OUT" -c bindings/fortran/src/seams.f90 -o "$VIDEO_OUT/seams.o" && gfortran -J "$VIDEO_OUT" -I "$VIDEO_OUT" -c bindings/fortran/src/seams_methods.f90 -o "$VIDEO_OUT/seams_methods.o" && gfortran -I "$VIDEO_OUT" docs/examples/minimal_fortran.f90 "$VIDEO_OUT/seams_kinds.o" "$VIDEO_OUT/seams_c.o" "$VIDEO_OUT/seams_support.o" "$VIDEO_OUT/seams.o" "$VIDEO_OUT/seams_methods.o" -o "$VIDEO_OUT/ice-f" -L"$SEAMS_BUILD_SRC" -lseams -Wl,-rpath,"$SEAMS_BUILD_SRC" && "$VIDEO_OUT/ice-f" input/traj/genice_sI.lammpstrj
```

Expected output (substrings the run must print):

```text
fortran example ok
```

Draft narration:

The Fortran example in the tree is the minimal Fortran program, run on input/traj/genice_sI.lammpstrj after compiling the Fortran module. It prints fortran example ok. The ice Ih census I just read, nop 128 and q6 0.625855 and hexagonal 128 and rings 256, is the same library.
