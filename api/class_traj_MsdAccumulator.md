<!-- https://docs.dseams.info/api/class_traj_MsdAccumulator.md -->

# class traj::MsdAccumulator {#class-traj-msdaccumulator}

## Overview {#overview}

Streaming multi-origin MSD accumulator with arbitrary-axis split. [More…](#details-classtraj-1-1-msd-accumulator)

## Detailed Documentation {#detailed-documentation}

Streaming multi-origin MSD accumulator with arbitrary-axis split.

When maxLag > 0, memory is bounded to O(nParticles * maxLag / originStride). When maxLag == 0, all lags up to framesSeen - 1 are retained.

### Methods {#methods}

Push one frame of coordinates.

If box9OrNull is non-null, coordinates are unwrapped incrementally across consecutive frames before accumulating displacements. If null, xyz is used directly as unwrapped coordinates.
