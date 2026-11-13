# Video overlays

Picture tracks are terminal sessions. Molecule stills that can sit over
those sessions are drawn with chemparseplot from the dumps in `input/`.

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

The PNG files are not committed. Odysee hosts the finished recordings.
`docs/videos/<slug>/embed.txt` is where the claim id goes once the
upload exists.
