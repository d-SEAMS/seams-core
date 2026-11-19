<!-- https://docs.dseams.info/api/group_capi_frame_io.md -->

# Frame reload and read-back {#frame-reload-and-read-back}

## Overview {#overview}

Overwrite coordinates in place and copy the stored frame out. [More…](#details-group-capi-frame-io)

## Detailed Documentation {#detailed-documentation}

Overwrite coordinates in place and copy the stored frame out.

There is no borrowed-pointer load. Every analysis reads the library’s per-atom records, so a view onto caller memory would be copied on each analysis call instead of once per frame. A trajectory loop calls [seams_system_set_atoms()](generated_targets.md#doxid-group-capi-frame-1ga483ad3783a7f233611be4b7d367a750f) once and [seams_system_set_positions()](generated_targets.md#doxid-group-capi-frame-io-1ga6d1b3f2819220c5e4081789d7612b572) for each later frame, which writes only the coordinates.
