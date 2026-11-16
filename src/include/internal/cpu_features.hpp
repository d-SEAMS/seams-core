// Copyright (c) d-SEAMS developers
#pragma once

// x86 feature checks. __builtin_cpu_supports needs __cpu_model from
// compiler-rt or libgcc, which some macOS toolchains do not link. CPUID
// can trap to the hypervisor, so the answers are read once.
// GCC and Clang keep the inline xgetbv. MSVC uses _xgetbv and __cpuidex.

#include <compiler_compat.hpp>

// clang-cl defines __x86_64__ and _MSC_VER. Clang's <cpuid.h> #defines
// __cpuid, and the MSVC STL later includes <intrin.h>, whose declaration
// void __cpuid(int[4], int) is then a macro invocation. Every _MSC_VER
// compiler, including clang-cl, uses the intrinsic instead.
#if defined(SEAMS_HOST_X86) && defined(_MSC_VER)
#include <intrin.h>
#elif defined(__x86_64__) || defined(__i386__)
#include <cpuid.h>
// Clang-cl defines both _MSC_VER and __clang__, so it takes this branch.
// cpuid.h's __cpuid macro then rejects the function declaration in the
// MSVC <intrin.h> that <chrono> includes.
#ifdef __cpuid
#undef __cpuid
#endif
#ifdef __cpuidex
#undef __cpuidex
#endif
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
#if defined(SEAMS_HOST_X86)
  unsigned a = 0, b = 0, c = 0, d = 0;
  unsigned xlo = 0;
#if defined(_MSC_VER)
  int regs[4] = {};
  __cpuidex(regs, 1, 0);
  a = static_cast<unsigned>(regs[0]);
  b = static_cast<unsigned>(regs[1]);
  c = static_cast<unsigned>(regs[2]);
  d = static_cast<unsigned>(regs[3]);
#else
  if (!__get_cpuid(1, &a, &b, &c, &d)) {
    return out;
  }
#endif
  const bool osxsave = (c & (1u << 27)) != 0;
  const bool avx = (c & (1u << 28)) != 0;
  const bool fma = (c & (1u << 12)) != 0;
  if (!osxsave || !avx) {
    return out;
  }
#if defined(_MSC_VER)
  xlo = static_cast<unsigned>(_xgetbv(0));
#else
  unsigned xhi = 0;
  __asm__ volatile("xgetbv" : "=a"(xlo), "=d"(xhi) : "c"(0));
  (void)xhi;
#endif
  // The OS saves the YMM state (bits 1, 2) and the ZMM state (bits 5-7).
  const bool ymm = (xlo & 0x6u) == 0x6u;
  const bool zmm = ymm && (xlo & 0xe0u) == 0xe0u;
#if defined(_MSC_VER)
  __cpuidex(regs, 7, 0);
  a = static_cast<unsigned>(regs[0]);
  b = static_cast<unsigned>(regs[1]);
  c = static_cast<unsigned>(regs[2]);
  d = static_cast<unsigned>(regs[3]);
#else
  if (!__get_cpuid_count(7, 0, &a, &b, &c, &d)) {
    return out;
  }
#endif
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

// false holds the kernels behind seams::ylm::avx512Kernels,
// symm::avx2Kernels and the AVX2 cutoff stash scan on their generic path;
// true leaves the choice to the CPU again. For tests.
void allowSimdKernels(bool allow);
bool simdKernelsAllowed();

} // namespace seams::cpu
