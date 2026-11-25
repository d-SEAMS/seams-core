//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <thread_team.hpp>

#include <algorithm>
#include <atomic>
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

#include <linkcell.h>

namespace seams {
namespace {

int &pinDepth() {
  thread_local int depth = 0;
  return depth;
}

// The plan and the saved masks below belong to the process. A call that
// arrives from a second host thread while one holds them would rebuild
// the plan under the first thread's team, so it runs unplaced instead.
std::atomic<bool> gHeld{false};

bool &holdsPlacement() {
  thread_local bool held = false;
  return held;
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
double gBusyAt = -1.0e9;
int gNarrow = 0;
bool gCrowded = false;
// The sharers and the slice of the last split. A family that only loses
// members keeps it: re-splitting while the others still run their
// current calls moves ranks, and two of them meet on one cpu.
int gHeldIds[68];
int gHeldN = 0;
int gHeldPlan[kMaxCpu];
int gHeldPlanN = 0;
int gHeldParts = 1;
bool gPinned = false;
bool gMarked = false;
constexpr double kCalm = 1.0;

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

bool readStatus(const char *path, char *buf, int cap);

// Threads of this process that are runnable now, the caller included. Our
// own workers can still be spinning from the last region.
int ownRunning() {
  DIR *tasks = opendir("/proc/self/task");
  if (tasks == nullptr) {
    return 1;
  }
  int running = 0;
  char path[512];
  char buf[1024];
  while (const dirent *ent = readdir(tasks)) {
    if (ent->d_name[0] < '0' || ent->d_name[0] > '9') {
      continue;
    }
    std::snprintf(path, sizeof(path), "/proc/self/task/%s/stat", ent->d_name);
    if (!readStatus(path, buf, sizeof(buf))) {
      continue;
    }
    const char *close = std::strrchr(buf, ')');
    if (close != nullptr && close[1] == ' ' && close[2] == 'R') {
      ++running;
    }
  }
  closedir(tasks);
  return running < 1 ? 1 : running;
}

// The two counts are separate snapshots. A wide team whose workers wait
// out time slices drops in and out of the run queue between them, so a
// look can read our own workers as somebody else's: beside two busy
// threads, ten runnable with six of them ours read as four others, and
// the team ran unplaced at full width for the next second. Our count is
// read on both sides of the machine count and the larger kept. A thread
// that is runnable only for a moment would narrow the team as well, so
// the lower of two looks decides.
int othersRunning() {
  int others = 1 << 30;
  for (int look = 0; look < 2; ++look) {
    const int before = ownRunning();
    const int running = procsRunning();
    const int after = ownRunning();
    others = std::min(others, std::max(0, running - std::max(before, after)));
  }
  return others;
}

// othersRunning, the lowest of three looks 2 ms apart. On an idle machine
// one look in a few read a thread that was runnable for a moment, or one
// of our own workers between two counts, and the team stayed narrowed to
// seven of eight cpus until the next look.
int othersSteady(int first) {
  int others = first;
  for (int again = 0; again < 2 && others > 0; ++again) {
    usleep(2000);
    others = std::min(others, othersRunning());
  }
  return others;
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
  int copies[64];
  int ncopies = 0;
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
    if (sibling && *nowners < 64) {
      owners[(*nowners)++] = pid;
    } else if (runningCopy && ncopies < 64) {
      copies[ncopies++] = pid;
    }
  }
  closedir(proc);
  // Siblings split the mask among themselves. A running copy outside the
  // family is runnable only part of the time, and counting it would move
  // every sibling's rank each time it came and went, so two siblings
  // would land on one cpu. Running copies count only for a process with
  // no siblings.
  if (*nowners == 0) {
    for (int i = 0; i < ncopies; ++i) {
      owners[(*nowners)++] = copies[i];
    }
  }
}

// The cpus the caller's team may use: every OpenMP place under a host
// binding, which has already pinned this thread to one of them, and the
// mask the first call saw otherwise.
const cpu_set_t &teamCpus() {
  static cpu_set_t cpus;
  static const bool known = [] {
    CPU_ZERO(&cpus);
#ifdef SEAMS_HAS_OPENMP
    if (omp_get_proc_bind() != omp_proc_bind_false) {
      for (int p = 0; p < omp_get_num_places(); ++p) {
        int ids[kMaxCpu];
        const int n = std::min(omp_get_place_num_procs(p), kMaxCpu);
        omp_get_place_proc_ids(p, ids);
        for (int i = 0; i < n; ++i) {
          if (ids[i] >= 0 && ids[i] < CPU_SETSIZE) {
            CPU_SET(ids[i], &cpus);
          }
        }
      }
    }
#endif
    if (CPU_COUNT(&cpus) == 0) {
      captureAllowed();
      for (int i = 0; i < gAllowedN; ++i) {
        CPU_SET(gAllowed[i], &cpus);
      }
    }
    return true;
  }();
  (void)known;
  return cpus;
}

// Seconds the cpus in `cpus` spent busy, from the per-cpu lines of
// /proc/stat. Time stolen by the hypervisor counts as busy.
bool busySeconds(const cpu_set_t &cpus, double &seconds) {
  static const double tick = static_cast<double>(sysconf(_SC_CLK_TCK));
  FILE *file = std::fopen("/proc/stat", "r");
  if (file == nullptr || !(tick > 0.0)) {
    if (file != nullptr) {
      std::fclose(file);
    }
    return false;
  }
  char line[512];
  long long busy = 0;
  int seen = 0;
  while (std::fgets(line, sizeof(line), file) != nullptr) {
    if (std::strncmp(line, "cpu", 3) != 0) {
      break;
    }
    if (line[3] < '0' || line[3] > '9') {
      continue;
    }
    char *at = line + 3;
    const long cpu = std::strtol(at, &at, 10);
    if (cpu < 0 || cpu >= CPU_SETSIZE || !CPU_ISSET(static_cast<int>(cpu), &cpus)) {
      continue;
    }
    long long field[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    for (long long &value : field) {
      value = std::strtoll(at, &at, 10);
    }
    // user nice system idle iowait irq softirq steal
    busy += field[0] + field[1] + field[2] + field[5] + field[6] + field[7];
    ++seen;
  }
  std::fclose(file);
  seconds = static_cast<double>(busy) / tick;
  return seen > 0;
}

double processSeconds() {
  timespec ts{};
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
  return static_cast<double>(ts.tv_sec) + static_cast<double>(ts.tv_nsec) * 1.0e-9;
}

void choosePlan(int width) {
  (void)width;
  gPlanN = 0;
  gNarrow = 0;
  if (gAllowedN < 2) {
    return;
  }
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
  std::sort(ids, ids + nids);
  if (gHeldN > 1 && nids <= gHeldN) {
    bool subset = true;
    for (int i = 0; i < nids && subset; ++i) {
      subset = std::binary_search(gHeldIds, gHeldIds + gHeldN, ids[i]);
    }
    if (subset) {
      for (int i = 0; i < gHeldPlanN; ++i) {
        gPlan[i] = gHeldPlan[i];
      }
      gPlanN = gHeldPlanN;
      gParts = gHeldParts;
      gDecided = true;
      gRescan = false;
      return;
    }
  }
  gHeldN = 0;
  if (nids < 2) {
    // Only a look that would narrow a full team is taken again, so a
    // crowded process does not sleep on the timed call that looks.
    int others = othersRunning();
    if (others > 0 && !gCrowded) {
      others = othersSteady(others);
    }
    gCrowded = others > 0;
    // Other runnable threads hold half the mask or more. A worker tied to
    // one of those cpus waits for that cpu's time slice at every barrier,
    // and so do the regions that place and restore the team: a 2 ms call
    // took 60 ms beside an 8-thread spinner. Unplaced, the scheduler moves
    // the workers to whichever cpu frees up and FairWidth narrows the
    // team. A busy neighbour that sleeps between bursts can read as idle
    // for a moment, so the team stays unplaced for kCalm before the next
    // look.
    if (2 * others >= gAllowedN) {
      gBusyAt = monoNow();
      gParts = 1;
      gRescan = true;
      return;
    }
    // A few busy threads. A full team then has more workers than free
    // cpus and every barrier waits out a time slice: q6 through the frame
    // took 65 ms beside two busy threads at eight workers and 3 ms at six.
    // The team takes the free cpus and stays unplaced, so the scheduler
    // keeps it off the busy ones, wherever those run.
    if (others > 0) {
      gNarrow = gAllowedN - others;
      gParts = 1;
      gDecided = true;
      gRescan = false;
      return;
    }
  }
  int parts = nids;
  if (parts < 1) {
    parts = 1;
  }
  if (parts > gAllowedN) {
    parts = gAllowedN;
  }
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
    gHeldPlan[i] = slice[i];
  }
  gPlanN = nslice;
  gParts = parts;
  gHeldPlanN = nslice;
  gHeldParts = parts;
  for (int i = 0; i < nids; ++i) {
    gHeldIds[i] = ids[i];
  }
  gHeldN = nids;
  gDecided = true;
  gRescan = false;
}

// A host that bound the team (OMP_PROC_BIND, OMP_PLACES, ...) has pinned
// the caller to its first place, and a Rayon pool started from there would
// put every worker on that one cpu. Hand linkcell the cpus the host gave
// the team instead, once per width. No cpu outside those places is used.
void shareHostPlaces() {
#ifdef SEAMS_HAS_OPENMP
  static int sharedWidth = 0;
  const int width = std::max(1, omp_get_max_threads());
  if (width < 2 || width == sharedWidth || omp_in_parallel()) {
    return;
  }
  sharedWidth = width;
  cpu_set_t places;
  CPU_ZERO(&places);
#pragma omp parallel num_threads(width)
  {
    cpu_set_t mine;
    CPU_ZERO(&mine);
    if (sched_getaffinity(0, sizeof(mine), &mine) == 0) {
#pragma omp critical(seams_host_places)
      CPU_OR(&places, &places, &mine);
    }
  }
  int cpus[kMaxCpu];
  int n = 0;
  for (int cpu = 0; cpu < CPU_SETSIZE && n < kMaxCpu; ++cpu) {
    if (CPU_ISSET(cpu, &places)) {
      cpus[n++] = cpu;
    }
  }
  if (n > 1) {
    lc_set_places(cpus, n);
  }
#endif
}

void applyPin(bool team) {
  if (pinOptOut()) {
    return;
  }
  if (hostOwnsAffinity()) {
    shareHostPlaces();
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
  if (monoNow() - gBusyAt < kCalm) {
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
  if (gNarrow > 0) {
#ifdef SEAMS_HAS_OPENMP
    gSavedOmp = width;
    if (team && gNarrow < width) {
      omp_set_num_threads(gNarrow);
      gRestoreOmp = true;
      if (askedWidth() == 0) {
        askedWidth() = width;
        gMarked = true;
      }
    }
#endif
    gPinned = false;
    gApplied = true;
    gBegan = monoNow();
    return;
  }
  if (gPlanN < 1) {
    return;
  }
  // A team narrower than the mask has free cpus to move to, and a worker
  // tied to one cpu waits out whatever else lands there. Beside a thread
  // that ran 5 ms in every 10 on one cpu, 4-thread q6 took 6.1 ms pinned
  // and 4.4 ms left to the scheduler, as on an idle machine.
  if (gParts == 1 && width < gPlanN) {
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
  if (team) {
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
  }
  lc_set_places(gPlan, gPlanN);
  gPinned = team;
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
  if (gPinned) {
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
  }
#ifdef SEAMS_HAS_OPENMP
  if (gRestoreOmp) {
    omp_set_num_threads(gSavedOmp);
    gRestoreOmp = false;
  }
#endif
  if (gMarked) {
    askedWidth() = 0;
    gMarked = false;
  }
  gApplied = false;
}
#endif

} // namespace

int &askedWidth() {
  thread_local int width = 0;
  return width;
}

FairState &fairState() {
  static FairState state;
  return state;
}

int &fairDepth() {
  thread_local int level = 0;
  return level;
}

bool placementHeld() {
#ifdef __linux__
  return holdsPlacement() && gApplied;
#else
  return false;
#endif
}

// Cpus' worth of load other processes put on teamCpus(), or -1 where it
// cannot be read or the caller has a single cpu. The first look counts
// their runnable threads (othersSteady); later looks average over at
// least 0.1 s.
static double othersLoad(int &ncpu) {
#ifdef __linux__
  static std::mutex lock;
  static bool looked = false;
  static double lastAt = 0.0;
  static double lastBusy = 0.0;
  static double lastOwn = 0.0;
  static double others = 0.0;
  const cpu_set_t &cpus = teamCpus();
  ncpu = CPU_COUNT(&cpus);
  if (ncpu < 2) {
    return -1.0;
  }
  const std::lock_guard<std::mutex> hold(lock);
  const double now = monoNow();
  if (!looked || now - lastAt >= 0.1) {
    double busy = 0.0;
    if (!busySeconds(cpus, busy)) {
      return -1.0;
    }
    const double own = processSeconds();
    if (looked) {
      others = std::max(0.0, (busy - lastBusy) - (own - lastOwn)) / (now - lastAt);
    } else {
      others = static_cast<double>(othersSteady(othersRunning()));
      looked = true;
    }
    lastAt = now;
    lastBusy = busy;
    lastOwn = own;
  }
  return others;
#else
  ncpu = 0;
  return -1.0;
#endif
}

int loadWidth(int base) {
  int ncpu = 0;
  if (base < 2) {
    return 0;
  }
  const double others = othersLoad(ncpu);
  if (others < 0.0) {
    return 0;
  }
  const int half = std::max(1, base / 2);
  if (2.0 * others >= static_cast<double>(ncpu)) {
    return half;
  }
  const int busyCpus = static_cast<int>(others + 0.5);
  return std::max(half, std::min(base, ncpu - busyCpus));
}

bool othersHoldHalf() {
  int ncpu = 0;
  const double others = othersLoad(ncpu);
  return others >= 0.0 && 2.0 * others >= static_cast<double>(ncpu);
}

PinSession::PinSession(bool team) : outer_(pinDepth()++ == 0) {
  if (outer_) {
    widenAffinity();
    bool idle = false;
    outer_ = gHeld.compare_exchange_strong(idle, true, std::memory_order_acquire);
  }
  if (outer_) {
    holdsPlacement() = true;
#ifdef __linux__
    applyPin(team);
#else
    (void)team;
#endif
  }
}

PinSession::~PinSession() {
  if (outer_) {
#ifdef __linux__
    restorePin();
#endif
    holdsPlacement() = false;
    gHeld.store(false, std::memory_order_release);
  }
  --pinDepth();
}

} // namespace seams
