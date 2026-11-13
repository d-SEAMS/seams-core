# d-SEAMS CHILL+ against mdapy CHILL+

Same three frames, published thresholds, cutoff 3.5 A, one thread.
mdapy's class names are aligned with the CHILL+ table: hexagonal ice,
cubic ice, interfacial ice, hydrate (clathrate), interfacial hydrate,
and other for liquid or unclassified. Thresholds were not moved.

## interface

d-SEAMS: cubic=1 hexagonal=141 water=309 interfacial=45 clathrate=2 interclathrate=14
mdapy: cubic_ice=1 hexagonal_ice=139 hydrate=6 interfacial_hydrate=12 interfacial_ice=47 other=307
per-atom disagreements: 20

| oxygen | d-SEAMS | mdapy | same shell |
|---:|---|---|---|
| 71 | hexagonal_ice | interfacial_ice | true |
| 108 | other | interfacial_hydrate | true |
| 122 | interfacial_ice | hexagonal_ice | true |
| 151 | hexagonal_ice | interfacial_ice | true |
| 174 | hexagonal_ice | interfacial_ice | true |
| 225 | interfacial_hydrate | other | true |
| 266 | interfacial_hydrate | hydrate | true |
| 281 | other | interfacial_hydrate | true |
| 308 | other | interfacial_hydrate | true |
| 320 | interfacial_hydrate | other | true |
| 353 | other | interfacial_hydrate | true |
| 365 | hydrate | other | true |
| 386 | interfacial_hydrate | other | true |
| 393 | interfacial_hydrate | hydrate | true |
| 413 | other | interfacial_hydrate | true |
| 435 | interfacial_hydrate | hydrate | true |
| 452 | other | hydrate | true |
| 492 | other | interfacial_hydrate | true |
| 499 | interfacial_hydrate | hydrate | true |
| 512 | interfacial_hydrate | other | true |

The 3.5 A neighbour set is the same for 20 of these 20 oxygens. The classes that swap are adjacent rows of the published table: liquid against interfacial hydrate, hexagonal against interfacial ice, and clathrate against hydrate. A correlation that sits on an edge falls one way in d-SEAMS and the other in mdapy. The edges were not moved.

## stacking

d-SEAMS: cubic=128 hexagonal=320 water=64
mdapy: cubic_ice=128 hexagonal_ice=320 other=64
per-atom disagreements: 0

Every oxygen has the same class.

## melted

d-SEAMS: hexagonal=369 water=114 interfacial=24 clathrate=1 interclathrate=4
mdapy: hexagonal_ice=365 interfacial_hydrate=3 interfacial_ice=24 other=120
per-atom disagreements: 17

| oxygen | d-SEAMS | mdapy | same shell |
|---:|---|---|---|
| 48 | interfacial_hydrate | other | true |
| 52 | hydrate | interfacial_hydrate | true |
| 53 | other | interfacial_hydrate | true |
| 81 | interfacial_ice | other | true |
| 171 | interfacial_hydrate | other | true |
| 186 | hexagonal_ice | interfacial_ice | true |
| 206 | interfacial_ice | hexagonal_ice | true |
| 212 | other | interfacial_hydrate | true |
| 272 | interfacial_ice | other | true |
| 278 | interfacial_ice | other | true |
| 304 | interfacial_hydrate | other | true |
| 314 | hexagonal_ice | interfacial_ice | true |
| 336 | hexagonal_ice | interfacial_ice | true |
| 345 | interfacial_hydrate | other | true |
| 431 | interfacial_ice | other | true |
| 458 | hexagonal_ice | interfacial_ice | true |
| 470 | hexagonal_ice | interfacial_ice | true |

The 3.5 A neighbour set is the same for 17 of these 17 oxygens. The classes that swap are adjacent rows of the published table: liquid against interfacial hydrate, hexagonal against interfacial ice, and clathrate against hydrate. A correlation that sits on an edge falls one way in d-SEAMS and the other in mdapy. The edges were not moved.
