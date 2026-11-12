#!/bin/sh
# Copyright (c) d-SEAMS developers
# Link an engine driver with the system C++ compiler.
# The module compiler's library path pulls a libgfortran and an MPI
# that do not resolve against this libc, so liblammps and libgromacs
# never finish linking.
exec env -u LIBRARY_PATH -u CPATH -u C_INCLUDE_PATH -u CPLUS_INCLUDE_PATH \
  PATH=/usr/bin:/bin /usr/bin/g++ "$@"
