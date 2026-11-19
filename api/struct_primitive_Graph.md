<!-- https://docs.dseams.info/api/struct_primitive_Graph.md -->

# struct primitive::Graph {#struct-primitive-graph}

## Overview {#overview}

This is a per-frame object, containing all the vertices for the particular frame, along with the vector of rings generated. [More…](#details-structprimitive-1-1-graph)

## Detailed Documentation {#detailed-documentation}

This is a per-frame object, containing all the vertices for the particular frame, along with the vector of rings generated.

Contains specifically the members:

- **pts** : Collection of vertices. The index of each should be the index according to the PointCloud.
- **rings** : A row-ordered vector of vectors for the rings generated, containing the indices (not IDs) of each member of the rings.

### Fields {#fields}

Collection of vertices. The index of each should be the same as that in pointCloud
