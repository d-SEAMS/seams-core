! Read the GenIce sI cell, then neighbours, q6, CHILL+, and rings.
program minimal_fortran
  use seams
  implicit none
  integer, parameter :: kmax = 64
  real(seams_real) :: flat(3 * kmax), lx, ly, lz, mean
  real(seams_real), allocatable :: xyz(:, :), q4(:), q6(:)
  integer(seams_int), allocatable :: labels(:), counts(:), columns(:), hist(:)
  integer(seams_int) :: n, i
  type(seams_frame) :: frame
  type(seams_neighbours) :: nb
  type(seams_rings) :: rings
  character(len=512) :: path
  integer :: argc

  argc = command_argument_count()
  if (argc < 1) stop 1
  call get_command_argument(1, path)
  call read_dump(trim(path), flat, n, lx, ly, lz)
  if (n /= 46) stop 2
  xyz = reshape(flat(1:3 * n), [3, n])
  call set_num_threads(1_seams_int)
  frame = new_frame()
  call frame%set_ortho(xyz, lx, ly, lz)
  nb = frame%neighbours()
  call nb%columns(columns)
  if (size(columns) < 150) stop 7
  call nb%close()
  call frame%q4q6(q4, q6)
  mean = 0.0_seams_real
  do i = 1, n
    mean = mean + q6(i)
  end do
  mean = mean / n
  call frame%chill_plus(labels, counts)
  if (counts(5) /= 46) stop 10
  rings = frame%rings(max_depth=6_seams_int)
  call rings%histogram(hist)
  if (hist(5) /= 48 .or. hist(6) /= 54) stop 13
  if (mean <= 0.5_seams_real) stop 14
  print *, "fortran example ok"
  call rings%close()
  call frame%close()

contains
  subroutine read_dump(filename, coords, natom, edge_x, edge_y, edge_z)
    character(len=*), intent(in) :: filename
    real(seams_real), intent(out) :: coords(*), edge_x, edge_y, edge_z
    integer(seams_int), intent(out) :: natom
    character(len=256) :: line
    integer :: unit, ios, id, typ, axis
    real(seams_real) :: lo, hi, x, y, z
    natom = 0
    open(newunit=unit, file=filename, status="old", action="read")
    do
      read(unit, "(A)", iostat=ios) line
      if (ios /= 0) exit
      if (index(line, "ITEM: BOX BOUNDS") == 1) then
        do axis = 1, 3
          read(unit, *) lo, hi
          if (axis == 1) edge_x = hi - lo
          if (axis == 2) edge_y = hi - lo
          if (axis == 3) edge_z = hi - lo
        end do
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
