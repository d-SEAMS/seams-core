<!-- https://docs.dseams.info/api/namespace_rdf.md -->

# namespace rdf {#namespace-rdf}

## Overview {#overview}

Three-dimensional partial RDF g_IJ(r). [More…](#details-namespacerdf)

## Detailed Documentation {#detailed-documentation}

Three-dimensional partial RDF g_IJ(r).

Pair distances use [gen::periodicDist](namespace_gen.md#doxid-namespacegen-1a4038dadde6016e7a738c2cb015344ee3), so a sheared dump cell uses the same minimum image as the neighbour lists. The ideal-gas denominator uses the triclinic volume from [nneigh::dumpVolume](namespace_nneigh.md#doxid-namespacenneigh-1ad17e972366865905372aa796e946a351) (det H), not the product of bound spans.

### Global Functions {#global-functions}

Histogram I-J pair distances and normalize to the dump-cell ideal gas.

Running CN: 4 pi rho_J int_0^r s^2 g_IJ(s) ds. rho_J = nJ / volume.

Running site-site CN: 4 pi rho_J integral_0^{r} s^2 g_IJ(s) ds per bin.

First minimum of g after the first maximum. Returns the bin index, or -1.

Site-site CN integrated to rMax. Same kernel as runningCN; no default rMax.

Site-site CN integrated to rMax. rhoJ is the number density of type J.
