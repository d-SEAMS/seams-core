Neighbours, CHILL, CHILL+, rings and cages can keep water oxygens by atom type or by molecule, and leave ions out of that graph.
`seams --mol water` uses the atom that appears once in a multi-atom molecule.
`seams --mol 1,5-12` keeps those molecule ids.
A type or molecule that is not in the frame is an error.
On a LAMMPS dump with multi-atom molecules and no `--type`, the same choice replaces the old preference for type 2, which labelled hydrogens or sodium.
`compute dseams` takes `type I` and `mol water`. A mixed frame that sets neither is an error.
CHILL+ then reports `ion-ice`, `ion-liquid` and `ion-front` when ions remain in the frame.
