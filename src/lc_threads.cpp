//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// Linkcell publishes lc_set_threads. A build of that library from before
// the symbol existed still has to link. These definitions are weak, so a
// strong copy in liblinkcell replaces them.

extern "C" __attribute__((weak)) void lc_set_threads(int) {}

extern "C" __attribute__((weak)) int lc_num_threads(void) { return 1; }
