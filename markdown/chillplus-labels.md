# Published CHILL+ and the OVITO variant

The default rule is Nguyen and Molinero, J. Phys. Chem. B 119, 9369
(2015). A bond is staggered when c_ij <= -0.8 and eclipsed when
-0.35 <= c_ij <= 0.25, both edges closed. Interfacial ice also reads
a neighbour's staggered count. The thresholds are not adjusted.

OVITO 3.17.2 and mdapy 2.0.0 share a different rule. On this machine
`ovito-chill+` matches mdapy on every atom. The published rule does
not, because of those edges and the neighbour test.

Ideal ice is unaffected: every bond sits far inside the staggered
side, so both rules label the whole crystal cubic.

Step 15677000 of the figshare nucleation trajectory
(doi:10.6084/m9.figshare.11448702, 4096 molecules, not a 13553-atom
frame):

| class | published | mdapy / OVITO |
| --- | ---: | ---: |
| cubic | 47 | 48 |
| hexagonal | 62 | 54 |
| water / other | 3589 | 3360 |
| interfacial ice | 240 | 478 |
| clathrate / hydrate | 20 | 21 |
| interclathrate / interfacial hydrate | 138 | 135 |

The gap is mostly water that OVITO calls interfacial.

The same split is in `discovery/v5-session.txt`. On step 15677000,
`ovito-chill+` matches mdapy and OVITO on all 4096 atoms and the
published rule differs on 397. On `input/traj/interface_f0.lammpstrj`
the match is 864/864 and the published rule differs on 46. q6 on both
frames agrees with mdapy and OVITO to 1e-14. The triclinic 110592
crystal is cubic on every atom under both rules.
