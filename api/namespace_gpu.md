<!-- https://docs.dseams.info/api/namespace_gpu.md -->

# namespace gpu {#namespace-gpu}

## Overview {#overview}

## Detailed Documentation {#detailed-documentation}

### Global Functions {#global-functions}

xyz is frame-major, atom, xyz.

box is frame-major three lengths when boxLow is null and nBox is 3. When boxLow is set or nBox is at least 6, box is one dump box (3 spans or 6 span+tilt) shared by every frame and the device walk uses lammpsBoxToLinkcell. rc is the cell-list cutoff; 5.5 matches the host k-NN candidate shell.
