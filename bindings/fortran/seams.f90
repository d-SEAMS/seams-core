! Copyright (c) d-SEAMS developers
!
! The modules are in src/. This one file compiles all of them, in the
! order they use each other, for a build that names a single source:
!   gfortran -c bindings/fortran/seams.f90
! fpm and meson compile src/ and must not also compile this file.
include "src/seams_kinds.f90"
include "src/seams_c.f90"
include "src/seams_support.f90"
include "src/seams.f90"
include "src/seams_methods.f90"
