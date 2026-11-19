# Recordings

Picture tracks are terminal sessions at 1920 by 1080. Finished recordings
are hosted on Odysee. `docs/videos/<slug>/embed.txt` holds one line,
`claim-name:full-claim-id`, after that upload exists. Until the line is
there, that slug renders no player and no Odysee URL. Claim ids are
never invented.

Captions are burned into the upload. The VTT file,
`docs/videos/<slug>/captions.vtt`, stays in the tree.

A later mix of narration audio uses ffmpeg with `-c:a aac -aac_coder fast`.

# Overlays

Molecule stills that can sit over those sessions are drawn with
chemparseplot from the dumps in `input/`.

From a checkout, with `libseams.so` already built:

```bash
pip install 'chemparseplot>=1.14' ovito ase
export SEAMS_LIBDIR=$PWD/build/src
python docs/scripts/render_molecules.py --out docs/videos
```

`pixi run render-molecules --out docs/videos` is the same command when
that task is on the path and the interpreter can import chemparseplot.
`--slug brine` draws one recording. The files land in
`<dir>/<slug>/overlay.png` (2400 by 1800) and `overlay-web.png`.

The PNG files are not committed.
