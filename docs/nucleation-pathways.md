# Nucleation pathways

One saved mW window and one short homogeneous quench, read with the
same four labels on every frame. The CHILL+ staggered and eclipsed
edges are the published edges. They are not adjusted.

## What is computed

The cutoff is 3.5 Å. Primitive rings are the Franzblau rings of that
graph through seven members. CHILL+ is applied twice: the published
Nguyen and Molinero rule (`chill+`, closed edges, four-neighbour
shell, interfacial neighbour test) and the OVITO rule
(`ovito-chill+`, open edges, q3 on the whole shell, no neighbour
test). q6-bar is the Lechner–Dellago average from
`chill::steinhardtQlCutoff` on the same cutoff list. A molecule is
counted as q6-bar ice when q6-bar is greater than 0.36. Diamond
symmetry labels come from `symm::diamondFromCsr`: 1 is cubic diamond
and 0 is hexagonal diamond.

The 0.36 cut sits between two reference frames that are not the
nucleation window. On `input/traj/mW_cubic.lammpstrj` every one of the
4096 molecules is published cubic, OVITO cubic and cubic diamond, and
the q6-bar median is 0.628539. On the 512-atom `input/traj/dump-l.lammpstrj`
frame the q6-bar median is 0.012006 and no molecule is above 0.36.
That cell is about 11.8 Å on a side, so its ring census is not a
liquid baseline. The cut is not moved to match the nucleation run.

A core molecule is cubic or hexagonal. An extended molecule is cubic,
hexagonal or interfacial. The nucleus size of a mask is the largest
connected component of the cutoff graph inside the mask. q6-bar does
not separate ice Ic from ice Ih. The stacking fraction of a q6-bar
cluster is the diamond label of its members, cubic divided by cubic
plus hexagonal. The same ratio is used for the CHILL+ clusters, with
the CHILL+ cubic and hexagonal counts.

## Time

`nucleation.lammpstrj` (figshare 10.6084/m9.figshare.11448702.v1, md5
`c28fc209a2854d2ba49e8f9b55accbe8`) has 148 frames of 4096 molecules.
The step index runs from 15677000 to 15824000 in steps of 1000. The
file does not store the integrator step. The mW deck in
`input/mw215/in.equilibrate` and the quench deck in
`input/nucleation/in.homogeneous` both use a 10 fs step. Under that
step the saved interval is 10 ps and the window is 1.47 ns. Frames
are the certain unit. Picoseconds below use 10 fs.

The window does not start from a liquid. At frame 1 the largest
published core cluster already has 59 molecules, and that series never
falls below 59.

## Half-rise

For a series x, the half-rise is the first frame at which x has gained
half of the increase from the first frame to the maximum. It is a
marker on this one growth window. It is not a committor and it is not
the free-energy critical size. The last uninterrupted climb to the
maximum sits within a few frames of the maximum, because the series
wiggles. That climb is printed in the summary as `onset_*` and is not
used as a nucleus size.

The six-ring signal is n6 / (n5 + n6). The ice signal is the ice-like
count divided by 4096.

## Does the ring share move first?

Inside this window, yes. The six-ring share goes from 0.552 at frame 1
(3462 six-rings, 2808 five-rings) to a maximum of 0.631 at frame 139.
Its half-rise is frame 94, 930 ps after the first frame. The ice-like
populations reach their own half-rise later:

| ice signal | half-rise lag | time |
| --- | ---: | ---: |
| published CHILL+, extended | 11 frames | 110 ps |
| published CHILL+, core | 14 frames | 140 ps |
| ovito-chill+, extended | 11 frames | 110 ps |
| ovito-chill+, core | 15 frames | 150 ps |
| q6-bar above 0.36 | 11 frames | 110 ps |
| diamond cubic or hexagonal | 8 frames | 80 ps |

A positive lag means the ring share moves first. The correlation of
the frame-to-frame changes is weak: 0.177 for the published extended
population at a lag of 6 frames, and between 0.18 and 0.26 for the
other masks. The lead is in the slow rise. There is no sharp
precursor spike. Using the raw six-ring count instead of the share
gives a 14-frame lead over the published extended population.

The fresh 97336-molecule liquid below has a six-ring share of 0.442
at 20 ps after the quench. The figshare window opens at 0.552, with a
CHILL+ cluster already present. The ring lead measured here is the
lead during growth, not a demonstration that rings appear before any
ice-like molecule.

## Largest cluster and stacking

The minimum of the largest cluster is the first frame for the
published core (59), the OVITO core (58) and the diamond label (2).
The extended and q6-bar series dip once just after the start: 146 at
frame 2, 166 at frame 3, and 113 at frame 4. From there the largest
cluster grows. On this window the published core nucleus is already at
least the size that goes on to grow, so the free-energy critical size
is not inside the saved frames. The smallest published core cluster in
the window, 59, is an upper bound only if that component is the one
that reaches 459. The run does not follow cluster identity.

The half-rise size is the largest cluster at the half-rise of that
same cluster series. Stacking is cubic / (cubic + hexagonal) inside
that component.

| definition | first size | half-rise frame | half-rise size | half-rise Ic fraction | maximum frame | maximum size | maximum Ic fraction |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| published core | 59 | 104 | 263 | 0.510 (134/129) | 143 | 459 | 0.566 (260/199) |
| published extended | 177 | 105 | 490 | 0.489 (135/141) | 139 | 752 | 0.571 (264/198) |
| ovito core | 58 | 103 | 256 | 0.516 (132/124) | 143 | 447 | 0.570 (255/192) |
| ovito extended | 200 | 105 | 528 | 0.498 (136/137) | 139 | 789 | 0.572 (261/195) |
| q6-bar | 131 | 105 | 429 | 0.562 (50/39) | 143 | 699 | 0.602 (80/53) |
| diamond | 2 | 102 | 91 | 0.538 (49/42) | 138 | 160 | 0.562 (90/70) |

The half-rise Ic fraction is about one half for every CHILL+ mask.
The maximum is slightly cubic, 0.57. The q6-bar fraction uses only
the members that diamond can label, 89 of 429 at the half-rise and
133 of 699 at the maximum, so it is a label on the ordered core of
the q6-bar cluster rather than on every member. Diamond itself never
labels more than 160 molecules. It sees the ordered interior, not the
CHILL+ cluster.

Counts at the ends of the window:

| frame | time | published cubic | published hexagonal | published interfacial | OVITO cubic | OVITO hexagonal | q6-bar ice | diamond cubic | diamond hexagonal | six-rings |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0 | 47 | 62 | 240 | 48 | 54 | 284 | 3 | 5 | 3462 |
| 148 | 1.47 ns | 269 | 217 | 353 | 268 | 210 | 749 | 103 | 70 | 3955 |

The population ratio cubic / (cubic + hexagonal) moves from 47/109 =
0.431 to 269/486 = 0.553 for the published rule, and from 48/102 =
0.471 to 268/478 = 0.561 for OVITO. The ice that accumulates in the
window is mixed and finishes slightly cubic. That agrees with the cage
census on this same file, where the molecules gained are mostly cubic
or shared by both cage types.

## How far the definition moves the answer

Shifts are definition B minus definition A at the half-rise of each
largest-cluster series. A negative size means B is smaller.

| A | B | size shift | frame shift | time shift |
| --- | --- | ---: | ---: | ---: |
| published core | OVITO core | -7 | -1 | -10 ps |
| published extended | OVITO extended | +38 | 0 | 0 |
| published core | q6-bar | +166 | +1 | +10 ps |
| published extended | q6-bar | -61 | 0 | 0 |
| published core | diamond | -172 | -2 | -20 ps |
| published extended | diamond | -399 | -3 | -30 ps |

Switching between the published rule and OVITO, while keeping the core
or the extended mask, moves the half-rise size by 7 or 38 molecules
and the half-rise time by at most 10 ps. Replacing the published core
with q6-bar adds 166 molecules. Replacing it with the diamond label
removes 172. The extended mask, which counts interfacial molecules as
part of the nucleus, is 227 molecules larger than the core at the
half-rise (490 against 263) and one frame later.

The two CHILL+ rules therefore do not move the induction marker of
this window by a meaningful amount. Counting the interface, or
switching to q6-bar or to diamond, does.

## Homogeneous quenches

`input/nucleation/regenerate.py` runs `input/nucleation/in.homogeneous`
with LAMMPS and the mW Stillinger–Weber file from `input/mw215`. The
integrator step is 10 fs. A melt at 350 K is quenched to 180 K. The
velocity seed is 20261111. Coordinate dumps stay outside the work
tree. LAMMPS was the 29 Aug 2024 tag, MANYBODY only, one MPI rank and
8 OpenMP threads.

### 97336 molecules

The cell is 46³. The melt is 2000 steps (20 ps) and the quench is 8000
steps (80 ps), dumped every 2000 steps. LAMMPS reported 6.875
steps/s on the quench. Five frames were saved, from the start of the
quench to 80 ps later. The per-frame table is
`docs/nucleation/mw-97336.tsv`.

The six-ring share and the count of published crystal molecules both
jump in the first 20 ps. The largest cluster does not.

| time after quench | six-ring share | published cubic + hexagonal | largest published core cluster | largest published extended cluster | q6-bar ice | largest q6-bar cluster | diamond labels |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.442 | 17 | 1 | 5 | 288 | 2 | 0 |
| 20 ps | 0.528 | 840 | 7 | 28 | 3570 | 20 | 5 |
| 40 ps | 0.549 | 1143 | 8 | 34 | 4719 | 29 | 6 |
| 60 ps | 0.561 | 1398 | 10 | 44 | 5756 | 32 | 8 |
| 80 ps | 0.569 | 1647 | 8 | 57 | 6499 | 35 | 9 |

There are five frames, so the half-rise lag is not reported: it
requires eight. At the 20 ps dump interval the ring share and the
dispersed crystal population rise together. They do not form one
nucleus. At 80 ps the published rule finds 1647 cubic or hexagonal
molecules and the largest connected core has 8. The extended cluster,
which adds interfacial molecules, reaches 57. OVITO finds 1448 cubic
or hexagonal molecules and 14868 interfacial molecules at that frame,
against 6438 interfacial for the published rule. The largest OVITO
extended cluster is 66, nine molecules above the published extended
cluster. The core clusters are the same size, 8. q6-bar labels 6499
molecules and its largest cluster has 35, none of which carry a
diamond label. Diamond labels 9 molecules in the whole cell.

This is a deep quench. It makes many small ice-like patches. It does
not produce a cluster on the scale of the 59-molecule core already
present at the first frame of the figshare window. The six-ring share
of this liquid, 0.442 at the start of the quench, is the baseline
under the figshare window, which opens at 0.552.

### 1000000 molecules

A 100³ cell was melted for 100 steps and quenched for 200 steps, with
one saved frame at step 200. That is 1 ps after the quench. The melt
had not settled: its last temperature was 245 K, not 350 K. LAMMPS
ran at 0.60 steps/s. The frame is a liquid check, not a pathway. The
published rule finds 110 cubic and 807 hexagonal molecules. The
largest core cluster has 3 molecules. Diamond labels one molecule.
q6-bar is above 0.36 for 3832 molecules, in a largest cluster of 8.
The six-ring share is 0.456. The row is
`docs/nucleation/mw-1000000.tsv`. A longer million-molecule quench is
the same script with a larger step count. The coordinates are not
committed.

## Files

The per-frame table is `docs/nucleation/figshare-4096.tsv`. The
summary lines are `docs/nucleation/figshare-4096.summary.tsv`. The
figure is `docs/nucleation/figshare-4096.png`. The 97336-molecule
quench is `docs/nucleation/mw-97336.tsv` and
`docs/nucleation/mw-97336.png`. The census program is
`nucleation-pathway`. The figures are drawn by
`scripts/nucleation_figures.py`.

```bash
nucleation-pathway nucleation.lammpstrj \
  --out docs/nucleation/figshare-4096.tsv \
  --summary docs/nucleation/figshare-4096.summary.tsv \
  --threads 8
python3 scripts/nucleation_figures.py docs/nucleation/figshare-4096.tsv \
  --out docs/nucleation/figshare-4096.png --ps-per-frame 10
```

The figshare dump is not in the repository. The first frame of it is
`input/traj/nucleation_f1.lammpstrj`.
