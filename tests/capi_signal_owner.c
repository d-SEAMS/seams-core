/* Who owns SIGABRT once libseams is loaded.
 * A plain build keeps the backward-cpp stack trace. A sanitized build leaves
 * the signal alone, so the sanitizer and libFuzzer see the crash.
 * Usage: capi_signal_owner plain|sanitized
 */
#define _POSIX_C_SOURCE 200809L

#include "seams.h"

#include <signal.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "usage: capi_signal_owner plain|sanitized\n");
    return 2;
  }
  seams_system *sys = seams_system_create();
  seams_system_destroy(sys);
  struct sigaction now;
  memset(&now, 0, sizeof now);
  if (sigaction(SIGABRT, NULL, &now) != 0) {
    perror("sigaction");
    return 1;
  }
  const int installed = (now.sa_flags & SA_SIGINFO) ? now.sa_sigaction != NULL
                                                    : now.sa_handler != SIG_DFL;
  if (strcmp(argv[1], "sanitized") == 0 && installed) {
    fprintf(stderr, "a sanitized libseams installed a SIGABRT handler\n");
    return 1;
  }
  if (strcmp(argv[1], "plain") == 0 && !installed) {
    fprintf(stderr, "libseams did not install the stack trace handler\n");
    return 1;
  }
  printf("SIGABRT %s\n", installed ? "handled" : "default");
  return 0;
}
