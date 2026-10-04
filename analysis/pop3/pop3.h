#ifndef ANALYSIS_POP3_H_
#define ANALYSIS_POP3_H_

/* POP CoE efficiencies for one OpenMP region, summed over every armed call.
 * Load balance is the mean useful time over the max useful time.
 * Communication efficiency is the max useful time over the region runtime.
 * Parallel efficiency is their product, the mean useful time over runtime.
 * IPC is retired instructions over unhalted cycles while that useful time
 * runs. A guest without a PMU leaves ipc_errno set and ipc at 0.
 */

#ifdef __cplusplus
extern "C" {
#endif

#define POP3_VIEW_THREADS 128

struct pop3_view {
  int threads;
  int regions;
  double runtime_s;
  double useful_mean_s;
  double useful_max_s;
  double load_balance;
  double communication_efficiency;
  double parallel_efficiency;
  int ipc_errno;
  unsigned long long instructions;
  unsigned long long cycles;
  double ipc;
  double useful_s[POP3_VIEW_THREADS];
};

/* Drop accumulated times. Call while disarmed. */
void pop3_reset(void);

/* Nonzero arms the GOMP_parallel wrapper. */
void pop3_arm(int on);

/* Open the hardware instruction and cycle events once. A failure is recorded
 * on the next snapshot and the timed region does not touch the counters. */
void pop3_probe(void);

/* Copy the armed accumulation. Call while disarmed. Returns 0. */
int pop3_snapshot(struct pop3_view *out);

#ifdef __cplusplus
}
#endif

#endif
