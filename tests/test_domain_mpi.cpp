// Run under mpiexec at any rank count: every rank's gathered rings are the
// whole frame's, and the reduced ring counts and per-atom types equal the
// threaded single-process summary.

#include "domain_frames.hpp"

#include <bop.hpp>
#include <domain.hpp>
#include <neighbours.hpp>

#include <mpi.h>

#include <cmath>
#include <cstdio>
#include <vector>

molSys::PointCloud<molSys::Point<double>, double> iceCube(int cells) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const double bond = 2.75;
  const double lattice = 4.0 * bond / std::sqrt(3.0);
  const double side = cells * lattice;
  const double basis[8][3] = {{0, 0, 0},         {0.5, 0.5, 0},      {0.5, 0, 0.5},
                              {0, 0.5, 0.5},     {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25},
                              {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75}};
  int id = 1;
  for (int i = 0; i < cells; i++) {
    for (int j = 0; j < cells; j++) {
      for (int k = 0; k < cells; k++) {
        for (const auto &b : basis) {
          molSys::Point<double> p;
          p.type = 1;
          p.atomID = id;
          p.molID = id;
          p.x = (i + b[0]) * lattice;
          p.y = (j + b[1]) * lattice;
          p.z = (k + b[2]) * lattice;
          cloud.pts.push_back(p);
          cloud.idIndexMap[id] = id - 1;
          id++;
        }
      }
    }
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  cloud.box = {side, side, side};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  chill::setSteinhardtAtomSplit(false);
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
      const auto want = seams::domain::summarize(sortedRowRings(c, cutoff, depth),
                                                 c.nop, depth);
      const auto red = seams::domain::reduceRings(c, cutoff, depth, MPI_COMM_WORLD);
      if (red.sizeCount != want.sizeCount) {
        std::fprintf(stderr, "rank %d: depth %d ring counts differ", rank, depth);
        for (std::size_t i = 0; i < red.sizeCount.size(); i++) {
          const long long w = i < want.sizeCount.size() ? want.sizeCount[i] : -1;
          std::fprintf(stderr, " [%zu]%lld/%lld", i, red.sizeCount[i], w);
        }
        std::fprintf(stderr, "\n");
        bad++;
      }
      if (rank == 0 && red.atomType != want.atomType) {
        std::fprintf(stderr, "rank 0: depth %d atom types differ\n", depth);
        bad++;
      }
      if (rank != 0 && !red.atomType.empty()) {
        bad++;
      }
    }
  }
  int nRanks = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
  {
    const double cutoff = 3.5;
    const auto base = iceCube(6);
    const auto rows = nneigh::neighListO(cutoff, base, 1);
    auto serial = base;
    chill::classifyBonds(serial, rows, chill::chillPlusRule(), false);
    chill::getIceTypePlusNoPrint(serial, rows, false, true);
    const auto qSerial = chill::steinhardtQl(base, rows, 6);
    auto got = base;
    chill::distributedChillPlus(got, cutoff, chill::chillPlusRule(), true, MPI_COMM_WORLD);
    const auto q = chill::distributedSteinhardt(base, cutoff, 6, MPI_COMM_WORLD);
    chill::BondClassifier ovito = chill::chillPlusRule();
    ovito.coordinationNumber = 0;
    ovito.closed = false;
    auto serialOvito = base;
    chill::classifyBonds(serialOvito, rows, ovito, false);
    chill::getIceTypePlusNoPrint(serialOvito, rows, false, false);
    auto gotOvito = base;
    chill::distributedChillPlus(gotOvito, cutoff, ovito, false, MPI_COMM_WORLD);
    const auto share = seams::domain::decompose(base, rank, nRanks, cutoff);
    for (std::size_t k = 0; k < share.local.size(); k++) {
      if (!share.owned[k]) {
        continue;
      }
      const int frame = share.local[k];
      if (got.pts[static_cast<std::size_t>(frame)].iceType !=
          serial.pts[static_cast<std::size_t>(frame)].iceType) {
        std::fprintf(stderr, "rank %d: chill label %d\n", rank, frame);
        bad++;
        break;
      }
      if (gotOvito.pts[static_cast<std::size_t>(frame)].iceType !=
          serialOvito.pts[static_cast<std::size_t>(frame)].iceType) {
        std::fprintf(stderr, "rank %d: ovito label %d\n", rank, frame);
        bad++;
        break;
      }
      const double dq = std::fabs(q.ql[static_cast<std::size_t>(frame)] -
                                  qSerial.ql[static_cast<std::size_t>(frame)]);
      const double db = std::fabs(q.qlBar[static_cast<std::size_t>(frame)] -
                                  qSerial.qlBar[static_cast<std::size_t>(frame)]);
      if (dq > 1e-15 || db > 1e-15) {
        std::fprintf(stderr, "rank %d: q6 %d %.3e %.3e\n", rank, frame, dq, db);
        bad++;
        break;
      }
    }
  }
  int anyBad = 0;
  MPI_Allreduce(&bad, &anyBad, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Finalize();
  return anyBad == 0 ? 0 : 1;
}
