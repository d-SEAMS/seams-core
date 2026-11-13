! Copyright (c) d-SEAMS developers
!
!> Numeric kinds for the Fortran API.
!!
!! Kinds come from `iso_fortran_env`. The C entry points hand
!! `double` and `int` buffers across `iso_c_binding`, so `seams_real`
!! must be `c_double` and `seams_int` must be `c_int`. A mismatch would
!! silently reinterpret memory. The assertions divide by zero in a
!! constant expression, which a conforming compiler rejects.
module seams_kinds
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr, c_null_ptr, &
    c_size_t, c_char, c_loc, c_f_pointer, c_associated, c_null_char
  implicit none
  private

  !> Working precision for coordinates, cutoffs, and order parameters.
  integer, parameter, public :: seams_real = real64
  !> Default integer kind for counts, labels, and the stat argument.
  integer, parameter, public :: seams_int = int32

  public :: c_double, c_int, c_ptr, c_null_ptr, c_size_t, c_char
  public :: c_loc, c_f_pointer, c_associated, c_null_char

  integer, parameter, private :: assert_real_is_c_double = 1/merge(1, 0, seams_real == c_double)
  integer, parameter, private :: assert_int_is_c_int = 1/merge(1, 0, seams_int == c_int)
end module seams_kinds
