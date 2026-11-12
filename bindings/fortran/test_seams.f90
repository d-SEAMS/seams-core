! Copyright (c) d-SEAMS developers
program test_seams
  use, intrinsic :: iso_c_binding
  use seams
  implicit none
  integer, parameter :: nrep = 2
  integer, parameter :: n = nrep * nrep * nrep * 8
  real(c_double) :: xyz(3 * n), q4(n), q6(n), length, a, fraction
  integer(c_int) :: labels(n), counts(9), code, i, cubic
  type(c_ptr) :: sys, nb, rings
  integer(c_int) :: hist(0:6)
  real(c_double) :: fcc(3, 4), basis(3, 2)
  integer :: ii, jj, kk, f, s, cursor

  fcc = reshape([0d0, 0d0, 0d0, 0.5d0, 0.5d0, 0d0, 0.5d0, 0d0, 0.5d0, 0d0, 0.5d0, 0.5d0], [3, 4])
  basis = reshape([0d0, 0d0, 0d0, 0.25d0, 0.25d0, 0.25d0], [3, 2])
  a = 2.75d0 * 4d0 / sqrt(3d0)
  cursor = 0
  do ii = 0, nrep - 1
    do jj = 0, nrep - 1
      do kk = 0, nrep - 1
        do f = 1, 4
          do s = 1, 2
            cursor = cursor + 1
            xyz(3 * cursor - 2) = (ii + fcc(1, f) + basis(1, s)) * a
            xyz(3 * cursor - 1) = (jj + fcc(2, f) + basis(2, s)) * a
            xyz(3 * cursor) = (kk + fcc(3, f) + basis(3, s)) * a
          end do
        end do
      end do
    end do
  end do
  length = nrep * a

  if (seams_set_num_threads(1) /= SEAMS_OK) stop 1
  if (seams_get_num_threads() /= 1) stop 2
  sys = seams_system_create()
  if (.not. c_associated(sys)) stop 3
  code = seams_system_set_ortho(sys, xyz, n, 0d0, 0d0, 0d0, length, length, length, &
                                c_null_ptr, c_null_ptr)
  if (code /= SEAMS_OK) stop 4
  code = seams_chill_plus(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, counts)
  if (code /= SEAMS_OK) stop 5
  cubic = 0
  do i = 1, n
    if (labels(i) == SEAMS_LABEL_CUBIC) cubic = cubic + 1
  end do
  if (cubic /= n) stop 6
  code = seams_chill_plus(sys, 3.5d0, 0, SEAMS_GRAPH_MUTUAL, SEAMS_CHILL_OVITO, labels, counts)
  if (code /= SEAMS_OK .or. counts(1) /= n) stop 7
  code = seams_q4q6(sys, 3.5d0, 0, q4, q6)
  if (code /= SEAMS_OK) stop 8
  if (sum(q6) / n <= 0.3d0) stop 9
  code = seams_neighbours_build(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, nb)
  if (code /= SEAMS_OK) stop 10
  call seams_neighbours_destroy(nb)
  code = seams_rings_find(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, 6, rings)
  if (code /= SEAMS_OK) stop 11
  code = seams_rings_histogram(rings, hist, 7)
  if (code /= SEAMS_OK .or. hist(6) <= 0) stop 12
  code = seams_rings_fraction(rings, 6, fraction)
  if (code /= SEAMS_OK .or. fraction <= 0.5d0) stop 13
  call seams_rings_destroy(rings)
  call seams_system_destroy(sys)
  print *, "fortran ok"
end program test_seams
