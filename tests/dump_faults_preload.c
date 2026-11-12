/* Copyright (c) d-SEAMS developers */
/* Test-side faults for pread and fstat. The library keeps calling libc. */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

static int pread_faults_left(void) {
  const char *n = getenv("DUMP_PREAD_EINTR");
  if (n == NULL || n[0] == '\0') {
    return 0;
  }
  return atoi(n);
}

ssize_t pread(int fd, void *buf, size_t count, off_t offset) {
  static ssize_t (*real_pread)(int, void *, size_t, off_t) = NULL;
  static int left = -1;
  if (real_pread == NULL) {
    real_pread = dlsym(RTLD_NEXT, "pread");
  }
  if (left < 0) {
    left = pread_faults_left();
  }
  if (getenv("DUMP_PREAD_EIO") != NULL) {
    errno = EIO;
    return -1;
  }
  if (left > 0) {
    left -= 1;
    errno = EINTR;
    return -1;
  }
  return real_pread(fd, buf, count, offset);
}

int fstat(int fd, struct stat *buf) {
  static int (*real_fstat)(int, struct stat *) = NULL;
  static int seen = 0;
  if (real_fstat == NULL) {
    real_fstat = dlsym(RTLD_NEXT, "fstat");
  }
  seen += 1;
  /* The parse records one successful stat, then the change check fails. */
  if (getenv("DUMP_FSTAT_EIO") != NULL && seen > 1) {
    errno = EIO;
    return -1;
  }
  return real_fstat(fd, buf);
}
