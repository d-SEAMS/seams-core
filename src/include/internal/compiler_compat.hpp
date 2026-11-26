// Copyright (c) d-SEAMS developers
#pragma once

// GNU target attributes and the inline xgetbv stay on GCC and Clang,
// including clang-cl. MSVC has no target attribute; those kernels take
// the scalar twin, which uses the same operation order.

#if defined(__GNUC__) || defined(__clang__)
#define SEAMS_TARGET(spec) __attribute__((target(spec)))
// Clang rejects the GCC optimize attribute. The kernels keep their
// default optimisation there; GCC still requests the named pass.
#if defined(__clang__)
#define SEAMS_OPTIMIZE(spec)
#else
#define SEAMS_OPTIMIZE(spec) __attribute__((optimize(spec)))
#endif
#define SEAMS_WEAK __attribute__((weak))
#define SEAMS_NOINLINE __attribute__((noinline))
#define SEAMS_ALWAYS_INLINE __attribute__((always_inline)) inline
#define SEAMS_COLD __attribute__((cold))
#define SEAMS_EXPORT_USED __attribute__((visibility("default"), used))
#define SEAMS_GNU_TARGET 1
#else
#define SEAMS_TARGET(spec)
#define SEAMS_OPTIMIZE(spec)
#define SEAMS_WEAK
#define SEAMS_NOINLINE
#define SEAMS_ALWAYS_INLINE __forceinline
#define SEAMS_COLD
#define SEAMS_EXPORT_USED
#endif

#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
#define SEAMS_HOST_X86 1
#endif
