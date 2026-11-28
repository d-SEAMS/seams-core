<!-- https://docs.dseams.info/user/tutorials/brine.md -->

# Brine and NaCl {#brine-and-nacl}

The recording follows an ice front through a short brine run, then this
864-water frame with nine ions.

`input/traj/nacl_front.lammpstrj` is the basal interface with nine ions
added. The dump has 873 atoms. Type 1 is the 864 waters. Type 3 is
chloride, five ions, and type 4 is sodium, four ions. The water
coordinates are the same as `input/traj/interface_f0.lammpstrj`.

`seams read` keeps the oxygen.

```bash
seams read input/traj/nacl_front.lammpstrj --type 1
```

```default
nop 864 frame 1 box 22.9883 21.6093 52.5093
```

## Water labels {#water-labels}

CHILL+ runs on the water, so an ion receives no ice label. The census
matches the ion-free frame.

```bash
seams chill-plus input/traj/nacl_front.lammpstrj --type 1 --cutoff 3.5 --layers --axis 2
```

```default
nop 864 clathrate 3 cubic 1 hexagonal 352 interClathrate 25 interfacial 49 water 434
phi-c 0.00283286 stack HHHHHHHH....C.
```

## Which side of the front {#which-side-of-the-front}

Slabs are 3.7 Å along $z$, measured from the box origin. A water is
crystal when its label is cubic, hexagonal, reCubic, or reHex. A slab
is ice when at least half of its water is crystal.

The chloride at $z = 12$ Å is in the slab that starts at 11.1 Å.
All 72 waters in that slab are crystal. The other eight ions are in
slabs whose crystal fraction is under one half.

| start (Å) | ions | side |
| --- | --- | --- |
| 11.1 | 1 chloride | ice |
| 33.3 | 1 chloride, 1 sodium | liquid |
| 37.0 | 1 chloride, 1 sodium | liquid |
| 44.4 | 1 chloride, 1 sodium | liquid |
| 48.1 | 1 chloride, 1 sodium | liquid |

The slab at 44.4 Å is on the liquid side of that cut and still holds
one crystal water, the same molecule counted on the ion-free interface.

![NaCl front. Water colored by published CHILL+, with sodium and chloride.](../../_images/overlay-web.png)

*NaCl front. Water colored by published CHILL+, with sodium and chloride. The chloride in the ice slab has no ice label.*
