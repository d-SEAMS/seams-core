/* An idle, resident OpenMP process for the omp-scaling contention rows.
 *
 * It holds a team of OMP_NUM_THREADS workers and enters an empty parallel
 * region every 50 microseconds, so between regions its workers sit in the
 * runtime's default wait (OMP_WAIT_POLICY unset), spinning before they
 * block. It does no useful work and runs until it is killed. */
#include <omp.h>
#include <stdio.h>
#include <time.h>

int main(void) {
  volatile long sink = 0;
  const struct timespec pause = {0, 50000};
  printf("spinner threads=%d\n", omp_get_max_threads());
  fflush(stdout);
  for (;;) {
#pragma omp parallel
    {
      sink += omp_get_thread_num() == 0;
    }
    nanosleep(&pause, NULL);
  }
  return 0;
}
