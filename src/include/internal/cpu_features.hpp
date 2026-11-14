// Copyright (c) d-SEAMS developers
#pragma once

// x86 feature checks. __builtin_cpu_supports needs __cpu_model from
// compiler-rt or libgcc, which some macOS toolchains do not link. CPUID
// can trap to the hypervisor, so the answers are read once.

#if defined(__x86_64__) || defined(__i386__)
#include <cpuid.h>
#endif

namespace seams::cpu {

struct Features {
  bool avx2 = false;
  bool fma = false;
  bool avx512f = false;
  bool avx512dq = false;
  bool avx512vnni = false;
};

inline Features detectFeatures() {
  Features out;
#if defined(__x86_64__) || defined(__i386__)
  unsigned a = 0, b = 0, c = 0, d = 0;
  if (!__get_cpuid(1, &a, &b, &c, &d)) {
    return out;
  }
  const bool osxsave = (c & (1u << 27)) != 0;
  const bool avx = (c & (1u << 28)) != 0;
  const bool fma = (c & (1u << 12)) != 0;
  if (!osxsave || !avx) {
    return out;
  }
  unsigned xlo = 0, xhi = 0;
  __asm__ volatile("xgetbv" : "=a"(xlo), "=d"(xhi) : "c"(0));
  (void)xhi;
  // The OS saves the YMM state (bits 1, 2) and the ZMM state (bits 5-7).
  const bool ymm = (xlo & 0x6u) == 0x6u;
  const bool zmm = ymm && (xlo & 0xe0u) == 0xe0u;
  if (!__get_cpuid_count(7, 0, &a, &b, &c, &d)) {
    return out;
  }
  out.fma = ymm && fma;
  out.avx2 = ymm && (b & (1u << 5)) != 0;
  out.avx512f = zmm && (b & (1u << 16)) != 0;
  out.avx512dq = out.avx512f && (b & (1u << 17)) != 0;
  out.avx512vnni = out.avx512f && (c & (1u << 11)) != 0;
#endif
  return out;
}

inline const Features &features() {
  static const Features cached = detectFeatures();
  return cached;
}

} // namespace seams::cpu
