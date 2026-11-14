! Copyright (c) d-SEAMS developers
! SPDX-License-Identifier: MIT
! fpm test: cubic ice built in memory, then ice Ih from the repository
! when input/ is reachable (DSEAMS_ROOT, the repository root, or two
! levels up from bindings/fortran).
program check_dseams
  use, intrinsic :: iso_c_binding
  use seams
  implicit none
  integer, parameter :: nrep = 2
  integer, parameter :: n = nrep * nrep * nrep * 8
  real(c_double) :: xyz(3 * n), q4(n), q6(n), fcc(3, 4), basis(3, 2), a, length, fraction
  integer(c_int) :: labels(n), counts(9), hist(0:6), code, n_offsets, nnz
  type(c_ptr) :: sys, nb, rings, offsets, cols
  integer :: ii, jj, kk, f, s, cursor
  character(len=:), allocatable :: ice

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

  if (.not. c_associated(seams_version_string())) error stop "version"
  if (seams_set_num_threads(1) /= SEAMS_OK) error stop "threads"
  sys = seams_system_create()
  if (.not. c_associated(sys)) error stop "create"
  code = seams_system_set_ortho(sys, xyz, n, 0d0, 0d0, 0d0, length, length, length, &
                                c_null_ptr, c_null_ptr)
  if (code /= SEAMS_OK) error stop "set_ortho"

  code = seams_chill_plus(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, counts)
  if (code /= SEAMS_OK) error stop "chill_plus"
  if (count(labels == SEAMS_LABEL_CUBIC) /= n) error stop "cubic count"
  code = seams_q4q6(sys, 3.5d0, 0, q4, q6)
  if (code /= SEAMS_OK .or. sum(q6) / n <= 0.3d0) error stop "q6"
  code = seams_neighbours_build(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, nb)
  if (code /= SEAMS_OK) error stop "neighbours"
  code = seams_neighbours_csr(nb, offsets, n_offsets, cols, nnz)
  if (code /= SEAMS_OK .or. nnz /= 4 * n) error stop "csr"
  call seams_neighbours_destroy(nb)
  code = seams_rings_find(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, 6, rings)
  if (code /= SEAMS_OK) error stop "rings"
  code = seams_rings_histogram(rings, hist, 7)
  if (code /= SEAMS_OK .or. hist(6) <= 0) error stop "histogram"
  code = seams_rings_fraction(rings, 6, fraction)
  if (code /= SEAMS_OK .or. fraction <= 0.5d0) error stop "fraction"
  call seams_rings_destroy(rings)
  print '(a,i0,a,f8.6,a,i0)', 'dseams cubic ', count(labels == SEAMS_LABEL_CUBIC), &
    ' q6 ', sum(q6) / n, ' six-rings ', hist(6)

  ice = fixture('input/traj/ice1h.lammpstrj')
  if (len(ice) > 0) then
    code = read_frame(sys, ice, 1, -1)
    if (code /= SEAMS_OK .or. seams_system_count(sys) /= 128) error stop "ice Ih read"
    code = seams_chill_plus(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, &
                            labels, counts)
    if (code /= SEAMS_OK .or. counts(2) /= 128) error stop "ice Ih hexagonal"
    print '(a,i0)', 'dseams ice Ih hexagonal ', counts(2)
  else
    print '(a)', 'dseams input/ not found; ice Ih check skipped'
  end if
  call seams_system_destroy(sys)
  print '(a)', 'dseams fortran ok'

contains

  function fixture(rel) result(path)
    character(len=*), intent(in) :: rel
    character(len=:), allocatable :: path
    character(len=4096) :: root
    integer :: nroot, status
    logical :: there
    call get_environment_variable('DSEAMS_ROOT', root, nroot, status)
    if (status == 0 .and. nroot > 0) then
      path = root(1:nroot) // '/' // rel
      inquire (file=path, exist=there)
      if (there) return
    end if
    path = rel
    inquire (file=path, exist=there)
    if (there) return
    path = '../../' // rel
    inquire (file=path, exist=there)
    if (there) return
    path = ''
  end function

end program check_dseams
