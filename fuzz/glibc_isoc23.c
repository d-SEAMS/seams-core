// Copyright (c) d-SEAMS developers
/* The system libclang_rt.fuzzer archive is built on the runner's glibc.
 * glibc 2.38 redirects strto* to __isoc23_strto* in its headers, so that
 * archive references those symbols. The pixi gcc links an older libc that
 * does not provide them. Forwarding to the libc strto* this toolchain has
 * covers the integers libFuzzer parses. The definitions are weak so a
 * libc that already exports the symbol keeps its own copy. */

#include <stdlib.h>

__attribute__((weak)) long __isoc23_strtol(const char *nptr, char **endptr, int base) {
  return strtol(nptr, endptr, base);
}

__attribute__((weak)) unsigned long __isoc23_strtoul(const char *nptr, char **endptr,
                                                     int base) {
  return strtoul(nptr, endptr, base);
}

__attribute__((weak)) long long __isoc23_strtoll(const char *nptr, char **endptr, int base) {
  return strtoll(nptr, endptr, base);
}

__attribute__((weak)) unsigned long long __isoc23_strtoull(const char *nptr, char **endptr,
                                                           int base) {
  return strtoull(nptr, endptr, base);
}
