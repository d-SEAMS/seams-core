#!/bin/sh
# Strong and weak scaling on one release build.
# One Snakemake job, so a point's threads are the only threads on the host.
set -eu
meson setup bbdir-eessi --reconfigure -Dbuildtype=release -Dwith_cuda_neighbours=disabled
ninja -C bbdir-eessi tests/bench_scale tests/test_scale_match
meson test -C bbdir-eessi scale_match
snakemake -j1 scaling
