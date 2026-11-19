<!-- https://docs.dseams.info/api/namespace_sph.md -->

# namespace sph {#namespace-sph}

## Overview {#overview}

Functions used for spherical harmonics. [More…](#details-namespacesph)

## Detailed Documentation {#detailed-documentation}

Functions used for spherical harmonics.

The recommended function for calculating the spherical-harmonics based complex vector is [sph::spheriHarmo](#doxid-namespacesph_1aa9220b76da7400c9d934fa08ceeeba7a). This function uses the [Boost](https://www.boost.org/) libraries.

### Changelog {#changelog}

- Amrita Goswami \[[amrita16thaug646@gmail.com](mailto:amrita16thaug646@gmail.com)\]; date modified: Sept 19, 2019
- Rohit Goswami \[[rog32@hi.is](mailto:rog32@hi.is)\]; date modified: Mar 20, 2021

#### Global Functions {#global-functions}

Fill `out` with Y_lm, m = -l .. +l. Reuses `out` capacity.

Lookup table for Q3.

Lookup table for Q3 (m=0 to m=6).

Lookup table for Q4.

Lookup table for Q4 (m=0 to m=8).

Lookup table for Q6.

Lookup table for Q8.

Lookup table for Q8 (m=0 to m=16).

Lookup table for Q6 (m=0 to m=12).
