! Read the GenIce sI cell, then neighbours, q6, CHILL+, and rings.
program minimal_fortran
  use, intrinsic :: iso_c_binding
  use seams
  implicit none
  integer, parameter :: kmax = 64
  real(c_double) :: xyz(3 * kmax), q4(kmax), q6(kmax), lx, ly, lz, mean
  integer(c_int) :: labels(kmax), counts(9), hist(0:12), n, code, i
  integer(c_int) :: n_offsets, nnz
  type(c_ptr) :: sys, nb, rings, offsets, cols
  character(len=512) :: path
  integer :: argc

  argc = command_argument_count()
  if (argc < 1) stop 1
  call get_command_argument(1, path)
  call read_dump(trim(path), xyz, n, lx, ly, lz)
  if (n /= 46) stop 2
  if (seams_set_num_threads(1) /= SEAMS_OK) stop 3
  sys = seams_system_create()
  if (seams_system_set_ortho(sys, xyz, n, 0.0d0, 0.0d0, 0.0d0, lx, ly, lz, &
      c_null_ptr, c_null_ptr) /= SEAMS_OK) stop 4
  if (seams_neighbours_build(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, nb) /= SEAMS_OK) stop 5
  code = seams_neighbours_csr(nb, offsets, n_offsets, cols, nnz)
  if (code /= SEAMS_OK .or. nnz < 150) stop 7
  call seams_neighbours_destroy(nb)
  if (seams_q4q6(sys, 3.5d0, 0, q4, q6) /= SEAMS_OK) stop 8
  mean = 0.0d0
  do i = 1, n
    mean = mean + q6(i)
  end do
  mean = mean / n
  if (seams_chill_plus(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, &
      labels, counts) /= SEAMS_OK) stop 9
  if (counts(5) /= 46) stop 10
  if (seams_rings_find(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, 6, rings) /= SEAMS_OK) stop 11
  if (seams_rings_histogram(rings, hist, 13) /= SEAMS_OK) stop 12
  if (hist(5) /= 48 .or. hist(6) /= 54) stop 13
  if (mean <= 0.5d0) stop 14
  print *, "fortran example ok"
  call seams_rings_destroy(rings)
  call seams_system_destroy(sys)

contains
  subroutine read_dump(filename, coords, natom, edge_x, edge_y, edge_z)
    character(len=*), intent(in) :: filename
    real(c_double), intent(out) :: coords(*), edge_x, edge_y, edge_z
    integer(c_int), intent(out) :: natom
    character(len=256) :: line
    integer :: unit, ios, id, typ
    real(c_double) :: x, y, z
    natom = 0
    open(newunit=unit, file=filename, status="old", action="read")
    do
      read(unit, "(A)", iostat=ios) line
      if (ios /= 0) exit
      if (index(line, "ITEM: BOX BOUNDS") == 1) then
        block
          character(len=256) :: bound
          real(c_double) :: xy, xz, yz, xmin, xmax, ymin, ymax
          real(c_double) :: xlo_b, xhi_b, ylo_b, yhi_b, zlo_b, zhi_b
          integer :: ios
          xy = 0d0
          xz = 0d0
          yz = 0d0
          read(unit, "(A)") bound
          read(bound, *, iostat=ios) xlo_b, xhi_b, xy
          if (ios /= 0) read(bound, *) xlo_b, xhi_b
          read(unit, "(A)") bound
          read(bound, *, iostat=ios) ylo_b, yhi_b, xz
          if (ios /= 0) read(bound, *) ylo_b, yhi_b
          read(unit, "(A)") bound
          read(bound, *, iostat=ios) zlo_b, zhi_b, yz
          if (ios /= 0) read(bound, *) zlo_b, zhi_b
          xmin = min(0d0, xy, xz, xy + xz)
          xmax = max(0d0, xy, xz, xy + xz)
          ymin = min(0d0, yz)
          ymax = max(0d0, yz)
          edge_x = (xhi_b - xlo_b) - xmax + xmin
          edge_y = (yhi_b - ylo_b) - ymax + ymin
          edge_z = zhi_b - zlo_b
        end block
      else if (index(line, "ITEM: ATOMS") == 1) then
        do
          read(unit, *, iostat=ios) id, typ, x, y, z
          if (ios /= 0) exit
          natom = natom + 1
          coords(3 * natom - 2) = x
          coords(3 * natom - 1) = y
          coords(3 * natom) = z
        end do
      end if
    end do
    close(unit)
  end subroutine read_dump
end program minimal_fortran
