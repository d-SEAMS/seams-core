<!-- https://docs.dseams.info/api/struct_chill_VoronoiWeights.md -->

# struct chill::VoronoiWeights {#struct-chill-voronoiweights}

## Overview {#overview}

The Voronoi facet neighbours of one particle and their area weights. [More…](#details-structchill-1-1-voronoi-weights)

## Detailed Documentation {#detailed-documentation}

The Voronoi facet neighbours of one particle and their area weights.

Neighbours are those sharing a facet with the particle’s Voronoi cell; weights are facet areas normalized to sum to one. Cutoff-based parameters change discontinuously when a neighbour crosses the cutoff, which is the shortcoming Mickel, Kapfer, Schr”oder-Turk and Mecke identified (J. Chem. Phys. 138, 044501, 2013): facet areas shrink to zero continuously as a neighbour leaves the first shell.

### Fields {#fields}

Cloud indices of facet-sharing neighbours.

Facet areas, normalized to sum to one.
