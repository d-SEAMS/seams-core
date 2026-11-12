// One-session timing of the cutoff CSR kernels and the language front.
// Prints medians only. No host name.

#include <seams_lang.h>

#include <bop.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

int coreCount() {
#ifdef SEAMS_HAS_OPENMP
  return std::max(1, omp_get_num_procs());
#else
  return 1;
#endif
}

void setThreads(int n) {
#ifdef SEAMS_HAS_OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(n);
#else
  (void)n;
#endif
}

molSys::PointCloud<molSys::Point<double>, double> ice(int nrep) {
  const double a = 2.75 * 4.0 / std::sqrt(3.0);
  const double basis[2][3] = {{0, 0, 0}, {0.25, 0.25, 0.25}};
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const double length = nrep * a;
  cloud.box = {length, length, length};
  cloud.boxLow = {0, 0, 0};
  cloud.pts.reserve(static_cast<std::size_t>(nrep) * nrep * nrep * 8);
  for (int i = 0; i < nrep; ++i) {
    for (int j = 0; j < nrep; ++j) {
      for (int k = 0; k < nrep; ++k) {
        for (const auto &shift : fcc) {
          for (const auto &site : basis) {
            molSys::Point<double> p;
            p.x = (i + shift[0] + site[0]) * a;
            p.y = (j + shift[1] + site[1]) * a;
            p.z = (k + shift[2] + site[2]) * a;
            p.type = 1;
            p.atomID = cloud.nop + 1;
            cloud.idIndexMap[p.atomID] = cloud.nop;
            cloud.pts.push_back(p);
            cloud.nop += 1;
          }
        }
      }
    }
  }
  return cloud;
}

double median(std::vector<double> samples) {
  std::sort(samples.begin(), samples.end());
  const std::size_t n = samples.size();
  if (n == 0) {
    return 0.0;
  }
  if ((n % 2) == 1) {
    return samples[n / 2];
  }
  return 0.5 * (samples[n / 2 - 1] + samples[n / 2]);
}

double timeCsr(molSys::PointCloud<molSys::Point<double>, double> cloud, int threads,
               int samples) {
  setThreads(threads);
  std::vector<int> offsets;
  std::vector<int> cols;
  nneigh::orthoFullCsr(cloud, nullptr, 3.5, offsets, cols, false);
  chill::chillPlusFromCsr(cloud, offsets, cols);
  std::vector<double> ms;
  ms.reserve(static_cast<std::size_t>(samples));
  for (int s = 0; s < samples; ++s) {
    const auto t0 = std::chrono::steady_clock::now();
    nneigh::orthoFullCsr(cloud, nullptr, 3.5, offsets, cols, false);
    chill::chillPlusFromCsr(cloud, offsets, cols);
    const auto t1 = std::chrono::steady_clock::now();
    ms.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
  }
  return median(ms);
}

double timeLang(const molSys::PointCloud<molSys::Point<double>, double> &cloud, int threads,
                int samples) {
  setThreads(threads);
  const int n = cloud.nop;
  std::vector<double> xyz(static_cast<std::size_t>(n) * 3);
  std::vector<int> types(static_cast<std::size_t>(n), 1);
  std::vector<int> ids(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    xyz[static_cast<std::size_t>(3 * i)] = cloud.pts[static_cast<std::size_t>(i)].x;
    xyz[static_cast<std::size_t>(3 * i + 1)] = cloud.pts[static_cast<std::size_t>(i)].y;
    xyz[static_cast<std::size_t>(3 * i + 2)] = cloud.pts[static_cast<std::size_t>(i)].z;
    ids[static_cast<std::size_t>(i)] = i + 1;
  }
  const double length = cloud.box[0];
  const double box[9] = {0, 0, 0, length, length, length, 0, 0, 0};
  seams_frame *frame = seams_frame_from_arrays(xyz.data(), n, box, types.data(), ids.data());
  std::vector<int> labels(static_cast<std::size_t>(n));
  int counts[9] = {};
  seams_frame_chill_plus(frame, 3.5, 0, "cutoff", labels.data(), counts);
  std::vector<double> ms;
  ms.reserve(static_cast<std::size_t>(samples));
  for (int s = 0; s < samples; ++s) {
    const auto t0 = std::chrono::steady_clock::now();
    seams_frame_chill_plus(frame, 3.5, 0, "cutoff", labels.data(), counts);
    const auto t1 = std::chrono::steady_clock::now();
    ms.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
  }
  seams_frame_free(frame);
  return median(ms);
}

} // namespace

int main() {
  const int nrep = 24;
  auto cloud = ice(nrep);
  const int samples = 8;
  std::printf("atoms %d cores %d samples %d highway 1\n", cloud.nop, coreCount(), samples);
  std::printf("csr_chill_1thread_ms %.3f\n", timeCsr(cloud, 1, samples));
  std::printf("csr_chill_8thread_ms %.3f\n", timeCsr(cloud, 8, samples));
  std::printf("lang_chill_1thread_ms %.3f\n", timeLang(cloud, 1, samples));
  std::printf("lang_chill_8thread_ms %.3f\n", timeLang(cloud, 8, samples));
  return 0;
}
