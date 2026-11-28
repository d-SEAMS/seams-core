<!-- https://docs.dseams.info/api/struct_primitive_Vertex.md -->

# struct primitive::Vertex {#struct-primitive-vertex}

## Overview {#overview}

This is a collection of elements, for each point, required for graph traversal. [More…](#details-structprimitive-1-1-vertex)

## Detailed Documentation {#detailed-documentation}

This is a collection of elements, for each point, required for graph traversal.

Contains specifically the members:

- **atomIndex** : This is the index according to the PointCloud.
- **neighListIndex** : A vector of indices (not IDs) of the neighboring vertices.
- **inGraph** : Bool qualifier, which is true by default. Setting it to false removes the vertex from the graph.

### Fields {#fields}

This is the index according to pointCloud.

Contains the INDICES (not the atomIDs) of the neighbouring vertices
