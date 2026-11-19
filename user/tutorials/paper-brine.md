<!-- https://docs.dseams.info/user/tutorials/paper-brine.md -->

# Ice growing into brine {#ice-growing-into-brine}

The front is the height where the crystal slabs stop. CHILL+ labels
the waters. Slabs of 3.7 Å along z are written as a stack: `H` for
hexagonal, `C` for cubic, `.` for a slab that is neither. The front
is the first slab that is not crystal, times 3.7 Å. Ions are not
waters, so they stay out of that count. You still want them in the
box, because the liquid they sit in is the brine the ice is growing
into.

Packmol places six waters, one sodium and one chloride above a fixed
GenIce ice Ih block. OpenMM then relaxes three stages. At each later
stage the nearest brine oxygen is moved onto the next lattice plane
and restrained. That move is the growth step. The relax is why the
ions and the leftover liquid can answer it, instead of staying where
packmol left them.

The short example gives

```text
brine frame 0 hexagonal 31 water 38 front 7.40
brine frame 1 hexagonal 32 water 36 front 7.40
brine frame 2 hexagonal 32 water 35 front 7.40
```

The front stays at 7.40 Å, two slab widths, so the added hexagonal
water has not filled the next 3.7 Å slab. The hexagonal count rises
from 31 to 32 and stays there. The water count falls from 38 to 35.

## Full run {#full-run}

The supplied brine frame is
`input/traj/nacl_front.lammpstrj`: 864 waters and 9 ions, one frame.
CHILL+ is scored on the waters.

```text
brine full input molecules 873 frames 1 file-frames 1 file input/traj/nacl_front.lammpstrj
brine full paper waters 864 ions 9 frames 1
brine full waters 864 ions 9 hexagonal 352 water 434
```

Hexagonal is 352 and water is 434, the same counts as the
ion-free interface frame.
