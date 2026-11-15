//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// Linkcell publishes lc_set_threads. A build of that library from before
// the symbol existed still has to link. These definitions are weak, so a
// strong copy in liblinkcell replaces them.

#include <thread_team.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>

#ifdef __linux__
#include <dirent.h>
#include <fcntl.h>
#include <unistd.h>
#endif

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

#ifdef SEAMS_HAS_LINKCELL
#include <linkcell.h>
#endif

extern "C" __attribute__((weak)) void lc_set_places(const int *cpus, int n) {
  (void)cpus;
  (void)n;
}

extern "C" __attribute__((weak)) void lc_set_threads(int n) {
#ifdef SEAMS_HAS_OPENMP
  if (n > 0) {
    omp_set_num_threads(n);
  }
#else
  (void)n;
#endif
}

extern "C" __attribute__((weak)) int lc_num_threads(void) {
#ifdef SEAMS_HAS_OPENMP
  return std::max(1, omp_get_max_threads());
#else
  return 1;
#endif
}

namespace seams {
namespace {

int &pinDepth() {
  thread_local int depth = 0;
  return depth;
}

bool envSet(const char *name) {
  const char *value = std::getenv(name);
  return value != nullptr && value[0] != '\0';
}

bool pinOptOut() {
  const char *flag = std::getenv("SEAMS_PIN");
  return flag != nullptr && flag[0] == '0' && flag[1] == '\0';
}

bool startsWithToken(const char *text, const char *token) {
  const std::size_t n = std::strlen(token);
  if (std::strncmp(text, token, n) != 0) {
    return false;
  }
  return text[n] == '\0' || text[n] == ',';
}

bool hostOwnsAffinity() {
  if (envSet("OMP_PROC_BIND")) {
    const char *bind = std::getenv("OMP_PROC_BIND");
    if (std::strcmp(bind, "false") != 0 && std::strcmp(bind, "FALSE") != 0) {
      return true;
    }
  }
  if (envSet("OMP_PLACES") || envSet("GOMP_CPU_AFFINITY")) {
    return true;
  }
  const char *kmp = std::getenv("KMP_AFFINITY");
  if (kmp == nullptr || kmp[0] == '\0') {
    return false;
  }
  return !startsWithToken(kmp, "none") && !startsWithToken(kmp, "disabled");
}

#ifdef __linux__
constexpr int kMaxCpu = 256;

int gAllowed[kMaxCpu];
int gAllowedN = 0;
bool gAllowedKnown = false;
int gParts = 1;
bool gDecided = false;
bool gRescan = false;
int gCalls = 0;
bool gApplied = false;
int gUse = 0;
int gSavedOmp = 1;
bool gRestoreOmp = false;
cpu_set_t gSaved[kMaxCpu];
int gPlan[kMaxCpu];
int gPlanN = 0;

double monoNow() {
  timespec ts{};
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return static_cast<double>(ts.tv_sec) + static_cast<double>(ts.tv_nsec) * 1.0e-9;
}

double gBegan = 0.0;

void captureAllowed() {
  if (gAllowedKnown) {
    return;
  }
  gAllowedKnown = true;
  cpu_set_t set;
  CPU_ZERO(&set);
  if (sched_getaffinity(0, sizeof(set), &set) != 0) {
    return;
  }
  for (int cpu = 0; cpu < CPU_SETSIZE && gAllowedN < kMaxCpu; ++cpu) {
    if (CPU_ISSET(cpu, &set)) {
      gAllowed[gAllowedN++] = cpu;
    }
  }
}

int procsRunning() {
  FILE *file = std::fopen("/proc/stat", "r");
  if (file == nullptr) {
    return 1;
  }
  char line[256];
  int running = 1;
  while (std::fgets(line, sizeof(line), file) != nullptr) {
    if (std::strncmp(line, "procs_running ", 14) == 0) {
      running = std::atoi(line + 14);
      break;
    }
  }
  std::fclose(file);
  return running < 1 ? 1 : running;
}

bool readStatus(const char *path, char *buf, int cap) {
  const int fd = open(path, O_RDONLY);
  if (fd < 0) {
    return false;
  }
  const ssize_t n = read(fd, buf, static_cast<std::size_t>(cap - 1));
  close(fd);
  if (n <= 0) {
    return false;
  }
  buf[n] = '\0';
  return true;
}

const char *field(const char *text, const char *key) {
  const std::size_t n = std::strlen(key);
  const char *line = text;
  while (line != nullptr && *line != '\0') {
    if (std::strncmp(line, key, n) == 0) {
      const char *value = line + n;
      while (*value == ' ' || *value == '\t') {
        ++value;
      }
      return value;
    }
    line = std::strchr(line, '\n');
    if (line != nullptr) {
      ++line;
    }
  }
  return nullptr;
}

void ourIdentity(char *nameOut, int nameCap, int *ppid) {
  nameOut[0] = '\0';
  *ppid = -1;
  char buf[4096];
  if (!readStatus("/proc/self/status", buf, sizeof(buf))) {
    return;
  }
  const char *name = field(buf, "Name:");
  int n = 0;
  if (name != nullptr) {
    while (name[n] != '\0' && name[n] != '\n' && n < nameCap - 1) {
      nameOut[n] = name[n];
      ++n;
    }
  }
  nameOut[n] = '\0';
  const char *parent = field(buf, "PPid:");
  if (parent != nullptr) {
    *ppid = std::atoi(parent);
  }
}

bool sameName(const char *text, const char *want) {
  const char *name = field(text, "Name:");
  if (name == nullptr || want[0] == '\0') {
    return false;
  }
  const std::size_t n = std::strlen(want);
  return std::strncmp(name, want, n) == 0 && (name[n] == '\n' || name[n] == '\0');
}

void scanOthers(int *owners, int *nowners) {
  *nowners = 0;
  DIR *proc = opendir("/proc");
  if (proc == nullptr) {
    return;
  }
  const int self = static_cast<int>(getpid());
  char mine[64];
  int parent = -1;
  ourIdentity(mine, sizeof(mine), &parent);
  char path[512];
  char buf[4096];
  while (const dirent *ent = readdir(proc)) {
    if (ent->d_name[0] < '1' || ent->d_name[0] > '9') {
      continue;
    }
    const int pid = std::atoi(ent->d_name);
    if (pid == self || pid <= 0) {
      continue;
    }
    std::snprintf(path, sizeof(path), "/proc/%s/status", ent->d_name);
    if (!readStatus(path, buf, sizeof(buf))) {
      continue;
    }
    const char *ppid = field(buf, "PPid:");
    if (ppid != nullptr && std::atoi(ppid) == 2) {
      continue;
    }
    const char *state = field(buf, "State:");
    const int ppidValue = ppid != nullptr ? std::atoi(ppid) : -1;
    // A sibling of this program shares the mask even between calls.
    // Another copy that is actually on a cpu shares it too. Desktop
    // helpers and one-shot tools do not take a slice.
    const bool sibling = parent > 0 && ppidValue == parent && sameName(buf, mine);
    const bool runningCopy = state != nullptr && state[0] == 'R' && sameName(buf, mine);
    if ((sibling || runningCopy) && *nowners < 64) {
      owners[(*nowners)++] = pid;
    }
  }
  closedir(proc);
}

void choosePlan(int width) {
  (void)width;
  gPlanN = 0;
  if (gAllowedN < 2) {
    return;
  }
  const int running = procsRunning();
  int owners[64];
  int nowners = 0;
  // One status file per process, and only when a plan is being built.
  // A lone process pays this once. Siblings share a parent, so they are
  // visible even while they sit between calls.
  scanOthers(owners, &nowners);
  int ids[68];
  int nids = 0;
  ids[nids++] = static_cast<int>(getpid());
  for (int i = 0; i < nowners; ++i) {
      bool seen = false;
      for (int j = 0; j < nids; ++j) {
        if (ids[j] == owners[i]) {
          seen = true;
          break;
        }
      }
      if (!seen && nids < 68) {
        ids[nids++] = owners[i];
      }
  }
  int parts = nids;
  if (parts < 2 && running > 1 && !gDecided) {
    parts = running < gAllowedN ? running : gAllowedN;
  }
  if (parts < 1) {
    parts = 1;
  }
  if (parts > gAllowedN) {
    parts = gAllowedN;
  }
  std::sort(ids, ids + nids);
  int rank = static_cast<int>(static_cast<unsigned>(getpid()) % static_cast<unsigned>(parts));
  for (int i = 0; i < nids; ++i) {
    if (ids[i] == static_cast<int>(getpid())) {
      rank = i % parts;
      break;
    }
  }
  const int share = std::max(1, gAllowedN / parts);
  int slice[kMaxCpu];
  int nslice = 0;
  const int begin = rank * share;
  // The slice is a function of who is sharing, not of who got there
  // first. Two processes that see the same partners take disjoint cpus.
  for (int i = 0; i < share && begin + i < gAllowedN; ++i) {
    slice[nslice++] = gAllowed[begin + i];
  }
  if (nslice < 1) {
    gParts = parts;
    gDecided = true;
    gRescan = false;
    return;
  }
  for (int i = 0; i < nslice; ++i) {
    gPlan[i] = slice[i];
  }
  gPlanN = nslice;
  gParts = parts;
  gDecided = true;
  gRescan = false;
}

void applyPin() {
  if (pinOptOut() || hostOwnsAffinity()) {
    return;
  }
  captureAllowed();
  if (gAllowedN < 2) {
    return;
  }
#ifdef SEAMS_HAS_OPENMP
  if (omp_in_parallel()) {
    return;
  }
  const int width = std::max(1, omp_get_max_threads());
#else
  const int width = 1;
#endif
  if (width < 2) {
    return;
  }
  ++gCalls;
  // A lone process does not scan /proc on the timed calls. A shared
  // machine is noticed on the first busy entry and again every few
  // calls, which is what pulls two full-mask teams off the same cpus.
  const bool refresh = !gDecided || gRescan || gParts > 1 || (gCalls % 8) == 0;
  if (refresh) {
    if ((gCalls % 8) == 0) {
      gRescan = true;
    }
    choosePlan(width);
  }
  if (gPlanN < 1) {
    return;
  }
  const int use = std::min(width, gPlanN);
#ifdef SEAMS_HAS_OPENMP
  gSavedOmp = std::max(1, omp_get_max_threads());
  if (use != gSavedOmp) {
    omp_set_num_threads(use);
    gRestoreOmp = true;
  }
#endif
  gUse = use;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel num_threads(gUse)
#endif
  {
#ifdef SEAMS_HAS_OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    if (tid < kMaxCpu && sched_getaffinity(0, sizeof(gSaved[tid]), &gSaved[tid]) == 0) {
      cpu_set_t one;
      CPU_ZERO(&one);
      CPU_SET(gPlan[tid % gPlanN], &one);
      sched_setaffinity(0, sizeof(one), &one);
    }
  }
  lc_set_places(gPlan, gPlanN);
  gApplied = true;
  gBegan = monoNow();
}

void restorePin() {
  if (!gApplied) {
    return;
  }
  const double elapsed = monoNow() - gBegan;
  if (elapsed > 8.0e-3) {
    gRescan = true;
  }
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel num_threads(gUse)
#endif
  {
#ifdef SEAMS_HAS_OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    if (tid < gUse) {
      sched_setaffinity(0, sizeof(gSaved[tid]), &gSaved[tid]);
    }
  }
#ifdef SEAMS_HAS_OPENMP
  if (gRestoreOmp) {
    omp_set_num_threads(gSavedOmp);
    gRestoreOmp = false;
  }
#endif
  gApplied = false;
}
#endif

} // namespace

bool placementHeld() {
#ifdef __linux__
  return gApplied && pinDepth() > 0;
#else
  return false;
#endif
}

PinSession::PinSession() : outer_(pinDepth()++ == 0) {
  if (outer_) {
    widenAffinity();
#ifdef __linux__
    applyPin();
#endif
  }
}

PinSession::~PinSession() {
  if (outer_) {
#ifdef __linux__
    restorePin();
#endif
  }
  --pinDepth();
}

} // namespace seams
