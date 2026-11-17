//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

// Pick your poison.
//
// On GNU/Linux, you have few choices to get the most out of your stack trace.
//
// By default you get:
//	- object filename
//	- function name
//
// In order to add:
//	- source filename
//	- line and column numbers
//	- source code snippet (assuming the file is accessible)

// Install one of the following library then uncomment one of the macro (or
// better, add the detection of the lib and the macro definition in your build
// system)

// - apt-get install libdw-dev ...
// - g++/clang++ -ldw ...
// #define BACKWARD_HAS_DW 1

// - apt-get install binutils-dev ...
// - g++/clang++ -lbfd ...
// #define BACKWARD_HAS_BFD 1

// The stack walk needs <execinfo.h>, which musl does not have; there the
// crash handler is left out.
#if !defined(__linux__) || __has_include(<execinfo.h>)
#include "backward.hpp"

// A sanitizer owns the fatal signals. These handlers would replace its
// report, and libFuzzer would not install its own or keep the input.
#if defined(__SANITIZE_ADDRESS__) || defined(__SANITIZE_THREAD__)
#define SEAMS_SANITIZED 1
#elif defined(__has_feature)
#if __has_feature(address_sanitizer) || __has_feature(thread_sanitizer) ||           \
    __has_feature(memory_sanitizer) || __has_feature(undefined_behavior_sanitizer)
#define SEAMS_SANITIZED 1
#endif
#endif

#ifndef SEAMS_SANITIZED
namespace backward {

// Never destroyed. The handlers stay installed until the process ends, and
// a signal raised after static destruction (a sanitizer abort at exit) would
// otherwise read the freed alternate stack.
backward::SignalHandling &sh = *new backward::SignalHandling;

} // namespace backward
#endif // SEAMS_SANITIZED
#endif
