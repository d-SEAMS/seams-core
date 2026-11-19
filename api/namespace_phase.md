<!-- https://docs.dseams.info/api/namespace_phase.md -->

# namespace phase {#namespace-phase}

## Overview {#overview}

## Detailed Documentation {#detailed-documentation}

### Global Variables {#global-variables}

Mass densities, g/cm^3. Ih is the 1 bar crystal. Ice VI is the dense tetrahedral crystal, not ice I.

Ice/LDA split: the midpoint of the Ih and ice VI number densities. Ih (about 0.031) is ice. Ice VI (about 0.044) is not.

### Global Functions {#global-functions}

Open hexagonal channel along z: six-rings stacked, not a closed 512.

64-bit key of donated H directions per oxygen (molID-shared hydrogens).

Mean squared displacement of hydrogenType between two frames.

Local number density (1/A^3) in a sphere of radius rcut.

Bulk mass density in g/cm^3. The volume is the triclinic cell volume, not the product of the bound spans. molarMass is g/mol.

Water, 18.015 g/mol.

HDA/MDA/LDA/ice from local number density (1/A^3). Ice first-shell density is four neighbours in the 3.5 A sphere (~0.022 A^-3). HDA is the dense shell (>= 0.070 A^-3). ldaMax splits LDA from MDA.

Three-argument form: LDA/MDA split at the midpoint of iceMax and mdaMin.

Number density in 1/A^3 from a mass density in g/cm^3.

Local-density windows, 1/A^3. Ice is below kIceNumberDensityMax, LDA below 0.040, MDA below 0.070, else HDA.
