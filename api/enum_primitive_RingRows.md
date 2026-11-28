<!-- https://docs.dseams.info/api/enum_primitive_RingRows.md -->

# enum primitive::RingRows {#enum-primitive-ringrows}

Which adjacency the ring search packs before it walks. Auto uses one sorted array per row: a single pass on one thread, a team on several. Vectors keeps a std::vector per atom. Packed forces the array even when the team has one thread. The ring set does not depend on the choice.
