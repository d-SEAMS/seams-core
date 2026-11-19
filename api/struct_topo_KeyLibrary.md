<!-- https://docs.dseams.info/api/struct_topo_KeyLibrary.md -->

# struct topo::KeyLibrary {#struct-topo-keylibrary}

## Overview {#overview}

A dictionary from local keys to labels: the keys of reference structures (a polymorph library), so any molecule whose rooted neighbourhood matches a reference gets that reference’s name. [More…](#details-structtopo-1-1-key-library)

## Detailed Documentation {#detailed-documentation}

A dictionary from local keys to labels: the keys of reference structures (a polymorph library), so any molecule whose rooted neighbourhood matches a reference gets that reference’s name.

Keys are only comparable when they come from the same method, hop count and coloring, which the header records.

### Fields {#fields}

“nauty” or “wl”

keys were computed with vertex colors

key -> label
