! Neighbours, rings, and CHILL+ from Fortran.
! Every call is a status code. The frame pointer is released on the
! way out, including the error paths.
program lang_tasks
  use, intrinsic :: iso_c_binding
  implicit none

  interface
    function seams_frame_read(path, frame, type_filter) bind(C, name="seams_frame_read")
      import :: c_ptr, c_char, c_int
      type(c_ptr) :: seams_frame_read
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int), value :: frame, type_filter
    end function

    subroutine seams_frame_free(frame) bind(C, name="seams_frame_free")
      import :: c_ptr
      type(c_ptr), value :: frame
    end subroutine

    function seams_frame_nop(frame) bind(C, name="seams_frame_nop")
      import :: c_ptr, c_int
      integer(c_int) :: seams_frame_nop
      type(c_ptr), value :: frame
    end function

    function seams_frame_neighbours(frame, cutoff, type_filter, graph) &
        bind(C, name="seams_frame_neighbours")
      import :: c_ptr, c_double, c_int, c_char
      integer(c_int) :: seams_frame_neighbours
      type(c_ptr), value :: frame
      real(c_double), value :: cutoff
      integer(c_int), value :: type_filter
      character(kind=c_char), intent(in) :: graph(*)
    end function

    function seams_frame_neighbour_pairs(frame) bind(C, name="seams_frame_neighbour_pairs")
      import :: c_ptr, c_int
      integer(c_int) :: seams_frame_neighbour_pairs
      type(c_ptr), value :: frame
    end function

    function seams_frame_rings(frame, cutoff, type_filter, graph, max_depth, counts, counts_len) &
        bind(C, name="seams_frame_rings")
      import :: c_ptr, c_double, c_int, c_char
      integer(c_int) :: seams_frame_rings
      type(c_ptr), value :: frame
      real(c_double), value :: cutoff
      integer(c_int), value :: type_filter
      character(kind=c_char), intent(in) :: graph(*)
      integer(c_int), value :: max_depth
      integer(c_int), intent(out) :: counts(*)
      integer(c_int), value :: counts_len
    end function

    function seams_frame_chill_plus(frame, cutoff, type_filter, graph, labels, counts) &
        bind(C, name="seams_frame_chill_plus")
      import :: c_ptr, c_double, c_int, c_char
      integer(c_int) :: seams_frame_chill_plus
      type(c_ptr), value :: frame
      real(c_double), value :: cutoff
      integer(c_int), value :: type_filter
      character(kind=c_char), intent(in) :: graph(*)
      integer(c_int), intent(out) :: labels(*)
      integer(c_int), intent(out) :: counts(*)
    end function

    function seams_last_error() bind(C, name="seams_last_error")
      import :: c_ptr
      type(c_ptr) :: seams_last_error
    end function
  end interface

  character(len=512) :: path
  character(len=256) :: line
  character(kind=c_char), target :: graph(7)
  type(c_ptr) :: frame
  integer(c_int) :: pairs, status, nop
  integer(c_int) :: rings(7), counts(9)
  integer(c_int), allocatable :: labels(:)
  integer :: i
  character(len=256) :: root

  call get_environment_variable("DSEAMS_ROOT", root)
  if (len_trim(root) == 0) then
    write(0, '(A)') "DSEAMS_ROOT is unset"
    stop 1
  end if
  path = trim(root) // "/input/traj/ice1h.lammpstrj" // c_null_char
  graph = [character(kind=c_char) :: "c", "u", "t", "o", "f", "f", c_null_char]
  frame = seams_frame_read(path, 1, -1)
  if (.not. c_associated(frame)) then
    call report()
    stop 1
  end if
  status = seams_frame_neighbours(frame, 3.5d0, -1, graph)
  if (status /= 0) call bail(frame)
  pairs = seams_frame_neighbour_pairs(frame)
  rings = 0
  status = seams_frame_rings(frame, 3.5d0, -1, graph, 6, rings, 7)
  if (status /= 0) call bail(frame)
  nop = seams_frame_nop(frame)
  allocate(labels(nop))
  counts = 0
  status = seams_frame_chill_plus(frame, 3.5d0, -1, graph, labels, counts)
  if (status /= 0) call bail(frame)
  write(line, '(A,1X,I0,1X,A,7(1X,I0),1X,A,9(1X,I0))') &
    "pairs", pairs, "rings", (rings(i), i = 1, 7), "chill", (counts(i), i = 1, 9)
  write(*, '(A)') trim(line)
  call seams_frame_free(frame)

contains

  subroutine report()
    type(c_ptr) :: err
    character(kind=c_char), pointer :: text(:)
    err = seams_last_error()
    if (c_associated(err)) then
      call c_f_pointer(err, text, [128])
      write(0, '(128A)') (text(i), i = 1, 128)
    else
      write(0, '(A)') "seams call failed"
    end if
  end subroutine

  subroutine bail(held)
    type(c_ptr), value :: held
    call report()
    call seams_frame_free(held)
    stop 1
  end subroutine
end program
