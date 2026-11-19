<!-- https://docs.dseams.info/user/tutorials/gallery.md -->

# Science gallery {#science-gallery}

These pictures show the interface, brine, ring, cage, and nucleation examples.
Each caption explains the colors and reports the corresponding molecular
census.

## Interface slab {#interface-slab}

![Interface slab](../../_images/gallery-interface.png)

The picture is the 864-molecule paper frame, colored by published CHILL+. Published CHILL+ reports hexagonal 352 interfacial 49 water 434 clathrate 3 interClathrate 25. The OVITO rule reports hexagonal 346 interfacial 57 water 437 clathrate 3 interClathrate 20. Hexagonal ice is orange and liquid water is sky blue. The interfacial shell is where the two rules part.

```default
interface full input molecules 864 frames 1 file-frames 1 file input/traj/interface_f0.lammpstrj
interface full paper molecules 864 frames 1
interface full published hexagonal 352 interfacial 49 water 434 clathrate 3 interClathrate 25
interface full ovito hexagonal 346 interfacial 57 water 437 clathrate 3 interClathrate 20
interface full phi-c 0.00283286
interface full stack HHHHHHHH....C.
```

## Brine frame {#brine-frame}

![Brine frame](../../_images/gallery-brine.png)

The picture is the shipped brine frame, 864 waters and 9 ions. Sodium is vermillion and chloride is blue. Waters 864 ions 9 hexagonal 352 water 434.

```default
brine full input molecules 873 frames 1 file-frames 1 file input/traj/nacl_front.lammpstrj
brine full paper waters 864 ions 9 frames 1
brine full waters 864 ions 9 hexagonal 352 water 434
```

## Ice Ih six-rings {#ice-ih-six-rings}

![Ice Ih six-rings](../../_images/gallery-rings.png)

The picture is frame 1 of the 4096-molecule mW dump. The drawn bonds are the cutoff neighbours, and in this frame they close six-rings. The histogram is 6 8192. That frame has two six-rings per molecule.

```default
rings full input molecules 4096 frames 1 file-frames 11 file input/traj/mW_cubic.lammpstrj
rings full paper molecules 4096 frames 1 rings 8192
rings full nop 4096 rings 8192
rings full histogram 6 8192
```

## Nucleation nucleus {#nucleation-nucleus}

![Nucleation nucleus](../../_images/gallery-nucleation.png)

Each molecule is colored by its own q6bar on frame 2 of the 148-frame paper trajectory. Dark blue is liquid-like and yellow is ice-like. frame 1 q6bar 0.26063. frame 148 q6bar 0.287256.

```default
nucleation full input molecules 4096 frames 148 file-frames 148 file docs/.cache/figshare-11448702/nucleation.lammpstrj
nucleation full paper molecules 4096 frames 148
nucleation full frame 1 q6bar 0.26063
nucleation full frame 148 q6bar 0.287256
nucleation full lines 148
```

## sI and sII cages {#si-and-sii-cages}

![sI and sII cages](../../_images/gallery-cages-sI.png)

![sI and sII cages](../../_images/gallery-cages-sII.png)

The first cluster is sI and the second cluster is sII. Bonds draw the water rings that close the cages. sI holds 16 of the small cage and 48 of the large cage. sII holds 64 of the small cage and 32 of the large cage.

```default
cages full input sI molecules 368 frames 1 file sI-full.lammpstrj
cages full sI 16 48
cages full input sII molecules 544 frames 1 file sII-full.lammpstrj
cages full sII 64 32
cages full input sH molecules 272 frames 1 file sH-full.lammpstrj
cages full sH 24 16 8
cages full paper unit sI molecules 46 cages 2 6
cages full paper unit sII molecules 136 cages 16 8
cages full paper unit sH molecules 68 cages 6 4 2
```

The molecular views below highlight the structures discussed in the
tutorials.

### Molecule stills {#molecule-stills}

These views highlight the ice, ring, prism, cage, brine, interface, and nucleation structures used throughout the tutorials.

![Ice Ih colored by published CHILL+.](../../_images/overlay-web3.png)

*Ice Ih colored by published CHILL+.*

![Ice Ih with hydrogen atoms and the hydrogen-bond network.](../../_images/overlay-web4.png)

*Ice Ih with hydrogen atoms and the hydrogen-bond network.*

![Ice Ic colored by published CHILL+.](../../_images/overlay-web5.png)

*Ice Ic colored by published CHILL+.*

![Five-rings and six-rings drawn as translucent polygons.](../../_images/overlay-web.png)

*Five-rings and six-rings drawn as translucent polygons.*

![Upright six-ring faces of ice Ih.](../../_images/overlay-web6.png)

*Upright six-ring faces of ice Ih.*

![GenIce sI with a methane molecule in each cage.](../../_images/overlay-web2.png)

*GenIce sI with a methane molecule in each cage.*

![GenIce sII with methane in the small cages and THF in the large cages.](../../_images/overlay-web7.png)

*GenIce sII with methane in the small cages and THF in the large cages.*

![GenIce sH with methane and THF guests.](../../_images/overlay-web8.png)

*GenIce sH with methane and THF guests.*

![THF guest and the surrounding water.](../../_images/overlay-web9.png)

*THF guest and the surrounding water.*

![NaCl front. The chloride in the ice slice is not an ice label. CHILL+ runs on the water alone.](../../_images/overlay-web1.png)

*NaCl front. The chloride in the ice slice is not an ice label. CHILL+ runs on the water alone.*

![Basal interface, one ice bilayer. Crystal bonds stop at the ice. An interClathrate oxygen with only liquid neighbours is drawn as water.](../../_images/overlay-web10.png)

*Basal interface, one ice bilayer. Crystal bonds stop at the ice. An interClathrate oxygen with only liquid neighbours is drawn as water.*

![Nucleation frame colored by published CHILL+.](../../_images/overlay-web11.png)

*Nucleation frame colored by published CHILL+.*

![GenIce FAU host.](../../_images/overlay-web12.png)

*GenIce FAU host.*

![sI oxygen cell drawn with the solvis backend.](../../_images/overlay-web13.png)

*sI oxygen cell drawn with the solvis backend.*
