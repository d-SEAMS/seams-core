! Copyright (c) d-SEAMS developers
! Times the Fortran module over the C API. Caller arrays are allocated once.
program binding_overhead
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env, only: real64, int64
  use seams
  implicit none
  character(len=4096) :: line, path, op, held
  integer :: threads, n, ios
  integer(c_int) :: code
  integer(int64) :: ticks, t0, t1
  real(real64) :: rate, sink
  type(c_ptr) :: sys, nb, rings
  integer(c_int), allocatable :: labels(:)
  real(c_double), allocatable :: ql(:), bar(:)
  integer(c_int) :: counts(9), hist(7)
  real(c_double) :: fraction
  logical :: loaded

  call system_clock(count_rate=ticks)
  rate = real(ticks, real64)
  sink = 0.0_real64
  loaded = .false.
  n = 0
  held = ""
  sys = seams_system_create()
  path = ""
  do
    read(*, "(A)", iostat=ios) line
    if (ios /= 0) exit
    if (line(1:4) == "quit") exit
    call split_command(line, threads, op, path)
    call retune(threads)
    if (op /= "read" .and. (held /= path .or. .not. loaded)) then
      code = read_frame(sys, trim(path), 1, -1)
      call check(code, sys)
      n = seams_system_count(sys)
      held = path
      loaded = .true.
      call allocate_out(n, labels, ql, bar)
    end if
    if (op == "read") then
      call system_clock(t0)
      code = read_frame(sys, trim(path), 1, -1)
      call system_clock(t1)
      call check(code, sys)
      n = seams_system_count(sys)
      held = path
      loaded = .true.
      call allocate_out(n, labels, ql, bar)
      sink = sink + real(n, real64)
      call report(t0, t1, rate)
    else
      if (op == "knn") then
        call system_clock(t0)
        code = seams_neighbours_knn(sys, 4, 3.5d0, 0, 1, nb)
        call check(code, sys)
        call seams_neighbours_destroy(nb)
        call system_clock(t1)
        sink = sink + 1.0_real64
        call report(t0, t1, rate)
      else if (op == "steinhardt") then
        call system_clock(t0)
        code = seams_steinhardt(sys, 3.5d0, 0, 6, ql, int(n, c_size_t), bar, int(n, c_size_t))
        call system_clock(t1)
        call check(code, sys)
        sink = sink + ql(1) + bar(1)
        call report(t0, t1, rate)
      else if (op == "chill") then
        call system_clock(t0)
        code = seams_chill(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, &
            labels, int(n, c_size_t), counts)
        call system_clock(t1)
        call check(code, sys)
        sink = sink + real(counts(2), real64)
        call report(t0, t1, rate)
      else if (op == "rings") then
        call system_clock(t0)
        code = seams_rings_find(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, 6, rings)
        call check(code, sys)
        code = seams_rings_histogram(rings, hist, 7)
        call check(code, rings)
        code = seams_rings_fraction(rings, 6, fraction)
        call check(code, rings)
        call seams_rings_destroy(rings)
        call system_clock(t1)
        sink = sink + real(hist(7), real64) + fraction
        call report(t0, t1, rate)
      else
        write(0, "(A)") "error unknown op"
        stop 1
      end if
    end if
  end do
  write(0, "(A,F0.6)") "sink ", sink
  call seams_system_destroy(sys)

contains

  subroutine split_command(line, nthreads, operation, filename)
    character(len=*), intent(in) :: line
    integer, intent(out) :: nthreads
    character(len=*), intent(out) :: operation, filename
    integer :: i, n, start
    n = len_trim(line)
    i = 1
    do while (i <= n .and. line(i:i) == " ")
      i = i + 1
    end do
    start = i
    do while (i <= n .and. line(i:i) /= " ")
      i = i + 1
    end do
    read(line(start:i - 1), *) nthreads
    do while (i <= n .and. line(i:i) == " ")
      i = i + 1
    end do
    start = i
    do while (i <= n .and. line(i:i) /= " ")
      i = i + 1
    end do
    operation = line(start:i - 1)
    do while (i <= n .and. line(i:i) == " ")
      i = i + 1
    end do
    filename = line(i:n)
  end subroutine

  subroutine retune(nthreads)
    integer, intent(in) :: nthreads
    character(len=16) :: text
    integer(c_int) :: status
    write(text, "(I0)") nthreads
    call setenv("OMP_NUM_THREADS", trim(text), 1, status)
    call setenv("PLUMED_NUM_THREADS", trim(text), 1, status)
    if (seams_set_num_threads(int(nthreads, c_int)) /= SEAMS_OK) then
      write(0, "(A)") "error threads"
      stop 1
    end if
  end subroutine

  subroutine check(status, handle)
    integer(c_int), intent(in) :: status
    type(c_ptr), intent(in) :: handle
    type(c_ptr) :: text
    character(kind=c_char), pointer :: raw(:)
    if (status == SEAMS_OK) return
    text = seams_handle_error(handle)
    if (c_associated(text)) then
      call c_f_pointer(text, raw, [256])
      write(0, "(A)") "error call"
    else
      write(0, "(A)") "error call"
    end if
    stop 1
  end subroutine

  subroutine allocate_out(count, labels_out, ql_out, bar_out)
    integer, intent(in) :: count
    integer(c_int), allocatable, intent(out) :: labels_out(:)
    real(c_double), allocatable, intent(out) :: ql_out(:), bar_out(:)
    allocate(labels_out(count), ql_out(count), bar_out(count))
  end subroutine

  subroutine report(t0, t1, clock_rate)
    integer(int64), intent(in) :: t0, t1
    real(real64), intent(in) :: clock_rate
    write(*, "(F0.9)") real(t1 - t0, real64) / clock_rate * 1.0d3
    flush(6)
  end subroutine

  subroutine setenv(name, value, overwrite, status)
    character(len=*), intent(in) :: name, value
    integer, intent(in) :: overwrite
    integer, intent(out) :: status
    interface
      function c_setenv(name, value, overwrite) bind(C, name="setenv")
        use iso_c_binding
        character(kind=c_char), intent(in) :: name(*), value(*)
        integer(c_int), value :: overwrite
        integer(c_int) :: c_setenv
      end function
    end interface
    character(kind=c_char), allocatable :: nbuf(:), vbuf(:)
    nbuf = c_chars(name)
    vbuf = c_chars(value)
    status = c_setenv(nbuf, vbuf, overwrite)
  end subroutine
end program binding_overhead
