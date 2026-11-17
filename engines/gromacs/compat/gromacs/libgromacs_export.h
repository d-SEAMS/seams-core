/* Copyright (c) d-SEAMS developers
 * CMake generates this header when it builds GROMACS and does not install
 * it. The legacy trajectory-analysis headers named by gromacs_legacy_include
 * need it; libgromacs exports these symbols with default visibility. */
#ifndef SEAMS_LIBGROMACS_EXPORT_H
#define SEAMS_LIBGROMACS_EXPORT_H
#define LIBGROMACS_EXPORT __attribute__((visibility("default")))
#define LIBGROMACS_NO_EXPORT __attribute__((visibility("hidden")))
#define LIBGROMACS_DEPRECATED __attribute__((__deprecated__))
#endif
