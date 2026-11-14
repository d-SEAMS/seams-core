! Copyright (c) d-SEAMS developers
! SPDX-License-Identifier: MIT
! Label one cubic ice cell (8 waters) with CHILL+ through the dseams package.
program downstream
  use, intrinsic :: iso_c_binding
  use seams
  implicit none
  real(c_double), parameter :: sites(3, 8) = reshape([ &
    0d0, 0d0, 0d0, 0.5d0, 0.5d0, 0d0, 0.5d0, 0d0, 0.5d0, 0d0, 0.5d0, 0.5d0, &
    0.25d0, 0.25d0, 0.25d0, 0.75d0, 0.75d0, 0.25d0, 0.75d0, 0.25d0, 0.75d0, &
    0.25d0, 0.75d0, 0.75d0], [3, 8])
  real(c_double) :: a, xyz(24)
  integer(c_int) :: labels(8), counts(9)
  type(c_ptr) :: sys

  a = 2.75d0 * 4d0 / sqrt(3d0)
  xyz = reshape(sites * a, [24])
  sys = seams_system_create()
  if (seams_system_set_ortho(sys, xyz, 8, 0d0, 0d0, 0d0, a, a, a, c_null_ptr, c_null_ptr) &
      /= SEAMS_OK) error stop "frame"
  if (seams_chill_plus(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, counts) &
      /= SEAMS_OK) error stop "chill+"
  call seams_system_destroy(sys)
  print '(a,i0,a)', 'downstream: ', count(labels == SEAMS_LABEL_CUBIC), ' of 8 waters are cubic ice'
end program downstream
