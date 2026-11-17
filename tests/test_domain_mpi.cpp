// Run under mpiexec at any rank count: every rank's gathered rings are the
// whole frame's, and the reduced ring counts and per-atom types equal the
// threaded single-process summary.

#include "domain_frames.hpp"
#include "slab_frames.hpp"

#include <bop.hpp>
#include <domain.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <string>
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
    for (const int depth : {4, 6, 8, 10}) {
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
    auto serialOvito = base;
    chill::classifyBonds(serialOvito, rows, ovito, false, false);
    chill::getIceTypePlusNoPrint(serialOvito, rows, false, false);
    auto gotOvito = base;
    chill::distributedChillPlus(gotOvito, cutoff, ovito, false, MPI_COMM_WORLD, false);
    const auto share = seams::domain::gatheredWork().share;
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
  auto trajOf = [](const std::string &name) {
    namespace fs = std::filesystem;
    fs::path dir = fs::current_path();
    for (int up = 0; up < 6; up++) {
      if (fs::exists(dir / "traj" / name)) {
        return (dir / "traj" / name).string();
      }
      if (fs::exists(dir / "input" / "traj" / name)) {
        return (dir / "input" / "traj" / name).string();
      }
      if (!dir.has_parent_path() || dir == dir.parent_path()) {
        break;
      }
      dir = dir.parent_path();
    }
    return std::string();
  };
  auto checkLabels = [&](molSys::PointCloud<molSys::Point<double>, double> base,
                         const char *name) {
    if (base.nop <= 0) {
      if (rank == 0) {
        std::fprintf(stderr, "%s: empty frame\n", name);
      }
      bad++;
      return;
    }
    const double cut = 3.5;
    const auto rows = nneigh::neighListO(cut, base, 1);
    auto serial = base;
    chill::classifyBonds(serial, rows, chill::chillPlusRule(), false);
    chill::getIceTypePlusNoPrint(serial, rows, false, true);
    auto got = base;
    chill::distributedChillPlus(got, cut, chill::chillPlusRule(), true, MPI_COMM_WORLD);
    const auto share = seams::domain::gatheredWork().share;
    for (std::size_t k = 0; k < share.local.size(); k++) {
      if (!share.owned[k]) {
        continue;
      }
      const int frame = share.local[k];
      if (got.pts[static_cast<std::size_t>(frame)].iceType !=
          serial.pts[static_cast<std::size_t>(frame)].iceType) {
        std::fprintf(stderr, "rank %d: %s label %d got %d serial %d\n", rank, name, frame,
                     static_cast<int>(got.pts[static_cast<std::size_t>(frame)].iceType),
                     static_cast<int>(serial.pts[static_cast<std::size_t>(frame)].iceType));
        bad++;
        break;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    chill::distributedChillPlus(got, cut, chill::chillPlusRule(), true, MPI_COMM_WORLD);
    double sample[6];
    for (int n = 0; n < 6; n++) {
      MPI_Barrier(MPI_COMM_WORLD);
      const auto t0 = std::chrono::steady_clock::now();
      chill::distributedChillPlus(got, cut, chill::chillPlusRule(), true, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      sample[n] = std::chrono::duration<double, std::milli>(
                      std::chrono::steady_clock::now() - t0)
                      .count();
    }
    if (rank == 0) {
      std::sort(sample, sample + 6);
      const double median = 0.5 * (sample[2] + sample[3]);
      const auto work = seams::domain::gatheredWork();
      std::printf("%s nop %d median %.3f ms halo %.3f\n", name, base.nop, median,
                  work.haloFraction);
    }
  };
  {
    const std::string liquid = trajOf("dump-l.lammpstrj");
    const std::string slab = trajOf("interface_f0.lammpstrj");
    const std::string cell = trajOf("genice_sI.lammpstrj");
    molSys::PointCloud<molSys::Point<double>, double> cloud;
    if (!liquid.empty()) {
      cloud = sinp::readLammpsTrjO(liquid, 1, cloud, 1);
      checkLabels(cloud, "liquid");
    } else {
      bad++;
    }
    if (!slab.empty()) {
      cloud = {};
      cloud = sinp::readLammpsTrjO(slab, 1, cloud, 1);
      checkLabels(cloud, "ice-water");
    } else {
      bad++;
    }
    if (!cell.empty()) {
      cloud = {};
      cloud = sinp::readLammpsTrjO(cell, 1, cloud, 1);
      const double lx = cloud.box.size() > 0 ? cloud.box[0] : 1.0;
      const double ly = cloud.box.size() > 1 ? cloud.box[1] : 1.0;
      const double lz = cloud.box.size() > 2 ? cloud.box[2] : 1.0;
      const int tile = 3;
      molSys::PointCloud<molSys::Point<double>, double> hydrate;
      hydrate.box = {lx * tile, ly * tile, lz * tile};
      hydrate.boxLow = cloud.boxLow;
      int id = 1;
      for (int ix = 0; ix < tile; ix++) {
        for (int iy = 0; iy < tile; iy++) {
          for (int iz = 0; iz < tile; iz++) {
            for (const auto &src : cloud.pts) {
              molSys::Point<double> p = src;
              p.atomID = id;
              p.molID = id;
              p.type = 1;
              p.x = src.x + ix * lx;
              p.y = src.y + iy * ly;
              p.z = src.z + iz * lz;
              const double mid = hydrate.box[2] * 0.5;
              if (p.z >= mid) {
                std::uint32_t s = static_cast<std::uint32_t>(id) * 1664525u + 1013904223u;
                p.x += (static_cast<int>(s % 1000) / 1000.0 - 0.5) * 3.5;
                s = s * 1664525u + 1013904223u;
                p.y += (static_cast<int>(s % 1000) / 1000.0 - 0.5) * 3.5;
                s = s * 1664525u + 1013904223u;
                p.z += (static_cast<int>(s % 1000) / 1000.0 - 0.5) * 2.0;
                auto wrap = [](double v, double L) {
                  v -= L * std::floor(v / L);
                  return v < 0.0 ? v + L : v;
                };
                p.x = wrap(p.x, hydrate.box[0]);
                p.y = wrap(p.y, hydrate.box[1]);
                p.z = wrap(p.z, hydrate.box[2]);
                if (p.z < mid) {
                  p.z += mid;
                  if (p.z >= hydrate.box[2]) {
                    p.z = hydrate.box[2] * 0.5 + 0.1;
                  }
                }
              }
              hydrate.pts.push_back(p);
              hydrate.idIndexMap[id] = id - 1;
              id++;
            }
          }
        }
      }
      hydrate.nop = static_cast<int>(hydrate.pts.size());
      checkLabels(hydrate, "hydrate-interface");
    } else {
      bad++;
    }
  }
  auto checkSlab = [&](molSys::PointCloud<molSys::Point<double>, double> base,
                       const char *name) {
    const double cut = 3.5;
    const auto rows = nneigh::neighListO(cut, base, 1);
    auto serial = base;
    chill::classifyBonds(serial, rows, chill::chillPlusRule(), false);
    chill::getIceTypePlusNoPrint(serial, rows, false, true);
    const auto qSerial = chill::steinhardtQl(base, rows, 6);
    auto got = base;
    chill::distributedChillPlus(got, cut, chill::chillPlusRule(), true, MPI_COMM_WORLD);
    const auto q = chill::distributedSteinhardt(base, cut, 6, MPI_COMM_WORLD);
    const auto share = seams::domain::gatheredWork().share;
    for (std::size_t k = 0; k < share.local.size(); k++) {
      if (!share.owned[k]) {
        continue;
      }
      const int frame = share.local[k];
      if (got.pts[static_cast<std::size_t>(frame)].iceType !=
          serial.pts[static_cast<std::size_t>(frame)].iceType) {
        std::fprintf(stderr, "rank %d: %s label %d got %d serial %d\n", rank, name, frame,
                     static_cast<int>(got.pts[static_cast<std::size_t>(frame)].iceType),
                     static_cast<int>(serial.pts[static_cast<std::size_t>(frame)].iceType));
        bad++;
        break;
      }
      if (static_cast<std::size_t>(frame) >= q.ql.size() ||
          static_cast<std::size_t>(frame) >= qSerial.ql.size()) {
        std::fprintf(stderr, "rank %d: %s q6 length\n", rank, name);
        bad++;
        break;
      }
      const double dq = std::fabs(q.ql[static_cast<std::size_t>(frame)] -
                                  qSerial.ql[static_cast<std::size_t>(frame)]);
      const double db = std::fabs(q.qlBar[static_cast<std::size_t>(frame)] -
                                  qSerial.qlBar[static_cast<std::size_t>(frame)]);
      if (dq > 1e-15 || db > 1e-15) {
        std::fprintf(stderr, "rank %d: %s q6 %d %.3e %.3e\n", rank, name, frame, dq, db);
        bad++;
        break;
      }
    }
    const auto want =
        seams::domain::summarize(sortedRowRings(base, cut, 6), base.nop, 6);
    for (int pass = 0; pass < 2; pass++) {
      const auto red = seams::domain::reduceRings(base, cut, 6, MPI_COMM_WORLD);
      if (red.sizeCount != want.sizeCount) {
        std::fprintf(stderr, "rank %d: %s rings pass %d\n", rank, name, pass);
        bad++;
      }
    }
    if (rank == 0) {
      std::printf("%s nop %d checked\n", name, base.nop);
    }
  };
  checkSlab(seams::frames::iceWaterSlab(24), "ice-water-110592");
  checkSlab(seams::frames::hydrateLiquidSlab(24, 24, 4), "hydrate-105984");
  checkSlab(seams::frames::iceWaterSlab(48), "ice-water-884736");
  checkSlab(seams::frames::hydrateLiquidSlab(40, 24, 20), "hydrate-883200");
  int anyBad = 0;
  MPI_Allreduce(&bad, &anyBad, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Finalize();
  return anyBad == 0 ? 0 : 1;
}
