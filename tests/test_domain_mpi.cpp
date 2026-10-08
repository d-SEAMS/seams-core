// Run under mpiexec: every rank's gathered rings are the whole frame's.

#include "domain_frames.hpp"

#include <domain.hpp>

#include <mpi.h>

#include <cstdio>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const double cutoff = 3.4;
  int bad = 0;
  for (const auto &c : {tiltedFrame(16000, 78.4, 0.0, 0.0, 0.0, 4),
                        tiltedFrame(16000, 78.4, 9.0, 5.0, -7.0, 5)}) {
    for (const int depth : {4, 6, 8}) {
      const auto got =
          seams::domain::gatherRings(c, cutoff, depth, MPI_COMM_WORLD);
      if (got.empty() || got != sortedRowRings(c, cutoff, depth)) {
        std::fprintf(stderr, "rank %d: depth %d rings differ\n", rank, depth);
        bad++;
      }
    }
  }
  int anyBad = 0;
  MPI_Allreduce(&bad, &anyBad, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Finalize();
  return anyBad == 0 ? 0 : 1;
}
