// Copyright (c) d-SEAMS developers
/* Who owns SIGABRT once libseams is loaded.
 * libseams installs no handler in either build: the host process (Julia, R,
 * Python, a sanitizer or libFuzzer) keeps its own. The backward-cpp stack
 * trace lives in the seams executables, which leave the signal to the
 * sanitizer in a sanitized build.
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
  if (strcmp(argv[1], "plain") != 0 && strcmp(argv[1], "sanitized") != 0) {
    fprintf(stderr, "usage: capi_signal_owner plain|sanitized\n");
    return 2;
  }
  if (installed) {
    fprintf(stderr, "a %s libseams installed a SIGABRT handler\n", argv[1]);
    return 1;
  }
  printf("SIGABRT %s\n", installed ? "handled" : "default");
  return 0;
}
