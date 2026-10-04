#define _GNU_SOURCE

#include "pop3.h"

#include <dlfcn.h>
#include <errno.h>
#include <linux/perf_event.h>
#include <omp.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <time.h>
#include <unistd.h>

/* libgomp in this build exports GOMP_parallel and no OMPT entry, so the
 * region wrapper is the hook. Useful time is CLOCK_THREAD_CPUTIME_ID across
 * the outlined function. That clock stops while a thread sleeps in a futex
 * and keeps running while libgomp spins, so a short active wait is counted
 * as useful. The master records CLOCK_MONOTONIC across the joined call.
 */

typedef void (*gomp_parallel_fn)(void (*fn)(void *), void *data,
                                 unsigned num_threads, unsigned flags);
typedef void (*outlined_fn)(void *);

static gomp_parallel_fn real_parallel;
static atomic_int armed;
static atomic_int team;
static atomic_int regions;
static atomic_int ipc_state; /* 0 unprobed, 1 open, -1 unavailable */
static atomic_int ipc_errno_store;
static atomic_ullong instructions;
static atomic_ullong cycles;
static double useful[POP3_VIEW_THREADS];
static double runtime_s;

static __thread int depth;
static __thread int fd_ins = -1;
static __thread int fd_cyc = -1;
static __thread int perf_ready;

struct frame {
  outlined_fn fn;
  void *data;
};

static long perf_open(unsigned long config) {
  struct perf_event_attr attr;
  memset(&attr, 0, sizeof(attr));
  attr.size = sizeof(attr);
  attr.type = PERF_TYPE_HARDWARE;
  attr.config = config;
  attr.disabled = 1;
  attr.exclude_kernel = 1;
  attr.exclude_hv = 1;
  return syscall(__NR_perf_event_open, &attr, 0, -1, -1, 0);
}

static void perf_begin(void) {
  if (atomic_load(&ipc_state) != 1)
    return;
  if (!perf_ready) {
    fd_ins = (int)perf_open(PERF_COUNT_HW_INSTRUCTIONS);
    fd_cyc = (int)perf_open(PERF_COUNT_HW_CPU_CYCLES);
    if (fd_ins < 0 || fd_cyc < 0) {
      atomic_store(&ipc_errno_store, errno);
      atomic_store(&ipc_state, -1);
      if (fd_ins >= 0)
        close(fd_ins);
      if (fd_cyc >= 0)
        close(fd_cyc);
      fd_ins = -1;
      fd_cyc = -1;
      return;
    }
    perf_ready = 1;
  }
  ioctl(fd_ins, PERF_EVENT_IOC_RESET, 0);
  ioctl(fd_cyc, PERF_EVENT_IOC_RESET, 0);
  ioctl(fd_ins, PERF_EVENT_IOC_ENABLE, 0);
  ioctl(fd_cyc, PERF_EVENT_IOC_ENABLE, 0);
}

static void perf_end(void) {
  unsigned long long ins = 0;
  unsigned long long cyc = 0;
  if (!perf_ready)
    return;
  ioctl(fd_ins, PERF_EVENT_IOC_DISABLE, 0);
  ioctl(fd_cyc, PERF_EVENT_IOC_DISABLE, 0);
  if (read(fd_ins, &ins, sizeof(ins)) == (ssize_t)sizeof(ins))
    atomic_fetch_add(&instructions, ins);
  if (read(fd_cyc, &cyc, sizeof(cyc)) == (ssize_t)sizeof(cyc))
    atomic_fetch_add(&cycles, cyc);
}

static void trampoline(void *arg) {
  struct frame *f = (struct frame *)arg;
  struct timespec a, b;
  int tid;
  int nt;
  depth++;
  perf_begin();
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &a);
  f->fn(f->data);
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &b);
  perf_end();
  depth--;
  tid = omp_get_thread_num();
  nt = omp_get_num_threads();
  if (nt > atomic_load(&team))
    atomic_store(&team, nt);
  if (tid >= 0 && tid < POP3_VIEW_THREADS) {
    useful[tid] += (double)(b.tv_sec - a.tv_sec) +
                   (double)(b.tv_nsec - a.tv_nsec) * 1e-9;
  }
}

void gomp_parallel_pop(void (*fn)(void *), void *data, unsigned num_threads,
                       unsigned flags);
__asm__(".symver gomp_parallel_pop,GOMP_parallel@@GOMP_4.0");

void gomp_parallel_pop(void (*fn)(void *), void *data, unsigned num_threads,
                       unsigned flags) {
  struct frame f;
  struct timespec a, b;
  if (!atomic_load(&armed) || depth > 0) {
    real_parallel(fn, data, num_threads, flags);
    return;
  }
  f.fn = fn;
  f.data = data;
  clock_gettime(CLOCK_MONOTONIC, &a);
  real_parallel(trampoline, &f, num_threads, flags);
  clock_gettime(CLOCK_MONOTONIC, &b);
  runtime_s += (double)(b.tv_sec - a.tv_sec) +
               (double)(b.tv_nsec - a.tv_nsec) * 1e-9;
  atomic_fetch_add(&regions, 1);
}

void pop3_reset(void) {
  int i;
  atomic_store(&team, 0);
  atomic_store(&regions, 0);
  atomic_store(&instructions, 0);
  atomic_store(&cycles, 0);
  runtime_s = 0.0;
  for (i = 0; i < POP3_VIEW_THREADS; i++)
    useful[i] = 0.0;
}

void pop3_arm(int on) { atomic_store(&armed, on ? 1 : 0); }

void pop3_probe(void) {
  int ins = (int)perf_open(PERF_COUNT_HW_INSTRUCTIONS);
  int cyc = (int)perf_open(PERF_COUNT_HW_CPU_CYCLES);
  if (ins < 0 || cyc < 0) {
    atomic_store(&ipc_errno_store, errno);
    atomic_store(&ipc_state, -1);
  } else {
    atomic_store(&ipc_errno_store, 0);
    atomic_store(&ipc_state, 1);
  }
  if (ins >= 0)
    close(ins);
  if (cyc >= 0)
    close(cyc);
}

int pop3_snapshot(struct pop3_view *out) {
  int n;
  int i;
  double sum = 0.0;
  double maxu = 0.0;
  unsigned long long ins;
  unsigned long long cyc;
  memset(out, 0, sizeof(*out));
  n = atomic_load(&team);
  if (n > POP3_VIEW_THREADS)
    n = POP3_VIEW_THREADS;
  out->threads = n;
  out->regions = atomic_load(&regions);
  out->runtime_s = runtime_s;
  for (i = 0; i < n; i++) {
    out->useful_s[i] = useful[i];
    sum += useful[i];
    if (useful[i] > maxu)
      maxu = useful[i];
  }
  if (n > 0)
    out->useful_mean_s = sum / (double)n;
  out->useful_max_s = maxu;
  if (maxu > 0.0)
    out->load_balance = out->useful_mean_s / maxu;
  if (runtime_s > 0.0) {
    out->communication_efficiency = maxu / runtime_s;
    out->parallel_efficiency = out->useful_mean_s / runtime_s;
  }
  ins = atomic_load(&instructions);
  cyc = atomic_load(&cycles);
  out->instructions = ins;
  out->cycles = cyc;
  out->ipc_errno = atomic_load(&ipc_errno_store);
  if (cyc > 0)
    out->ipc = (double)ins / (double)cyc;
  return 0;
}

__attribute__((constructor)) static void pop3_init(void) {
  real_parallel = (gomp_parallel_fn)dlsym(RTLD_NEXT, "GOMP_parallel");
  if (!real_parallel) {
    fprintf(stderr, "pop3: GOMP_parallel is not in the next object\n");
    abort();
  }
}
