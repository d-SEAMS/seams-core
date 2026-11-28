<!-- https://docs.dseams.info/api/namespace_clump.md -->

# namespace clump {#namespace-clump}

## Overview {#overview}

Clustering functions. This namespace contains functions that are used in the for clustering ice-like particles. [More…](#details-namespaceclump)

## Detailed Documentation {#detailed-documentation}

Clustering functions. This namespace contains functions that are used in the for clustering ice-like particles.

Particles are clustered according to Stoddard’s algorithm [(Stoddard J. Comp. Phys., 27, 291, 1977)](https://www.sciencedirect.com/science/article/pii/0021999178900116).

The largest ice cluster determines the largest cluster of ice-like particles, using a “linked list” created by the clustering algorithm. The determination of ice-like molecules can be done using the $q_6$ parameter, CHILL or CHILL+ parameters.

The largest ice cluster thus determined can be re-centered for ease of visualization.

### Changelog {#changelog}

- Amrita Goswami \[[amrita16thaug646@gmail.com](mailto:amrita16thaug646@gmail.com)\]; date modified: Nov 13, 2019
- Rohit Goswami \[[rog32@hi.is](mailto:rog32@hi.is)\]; date modified: Mar 20, 2021

#### Global Functions {#global-functions}

Finds the largest ice cluster.

Get the linked list of a cluster, given by iceCloud, for a single cluster. Required for cluster re-centering

Does the cluster analysis of ice particles in the system. Returns a pointCloud of the largest ice cluster. The neighbour list returned is BY INDEX of the largest ice cluster pointCloud.

Recenters the coordinates of a pointCloud.
