! Copyright (c) d-SEAMS developers
!
!> Shared handle ownership and the stat/errmsg channel.
!!
!! Not a public API. `seams` is the module programs use. A handle is a
!! reference-counted box so a derived-type assignment shares the C
!! pointer and the last finalizer destroys it. Copies are not
!! thread-safe: one handle is used by one thread at a time, which is
!! the rule in `seams.h`.
module seams_support
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, c_char, &
    c_null_char, c_size_t, c_int, c_f_pointer
  use seams_c, only: seams_handle_error, seams_last_error, seams_system_destroy, &
    seams_neighbours_destroy, seams_rings_destroy, seams_cages_destroy, &
    seams_fingerprint_destroy, seams_brine_destroy, SEAMS_OK
  implicit none
  private

  public :: handle_box, handle_acquire, handle_release, handle_share, handle_get
  public :: c_to_f, f_to_c, check_status, copy_key

  integer, parameter, public :: kind_system = 1
  integer, parameter, public :: kind_neighbours = 2
  integer, parameter, public :: kind_rings = 3
  integer, parameter, public :: kind_cages = 4
  integer, parameter, public :: kind_topo = 5
  integer, parameter, public :: kind_brine = 6

  type, public :: handle_box
    type(c_ptr) :: raw = c_null_ptr
    integer :: refcount = 0
    integer :: kind = 0
  end type handle_box

  interface
    function strlen(s) bind(C, name="strlen")
      import :: c_ptr, c_size_t
      type(c_ptr), value :: s
      integer(c_size_t) :: strlen
    end function
  end interface

contains

  !> Allocate a box that owns `raw`. The first owner holds refcount 1.
  subroutine handle_acquire(box, raw, kind)
    type(handle_box), pointer, intent(out) :: box
    type(c_ptr), intent(in) :: raw
    integer, intent(in) :: kind
    allocate(box)
    box%raw = raw
    box%refcount = 1
    box%kind = kind
  end subroutine handle_acquire

  !> Drop one owner. The last owner destroys the C handle.
  subroutine handle_release(box)
    type(handle_box), pointer, intent(inout) :: box
    if (.not. associated(box)) return
    box%refcount = box%refcount - 1
    if (box%refcount > 0) then
      box => null()
      return
    end if
    if (c_associated(box%raw)) then
      select case (box%kind)
      case (kind_system)
        call seams_system_destroy(box%raw)
      case (kind_neighbours)
        call seams_neighbours_destroy(box%raw)
      case (kind_rings)
        call seams_rings_destroy(box%raw)
      case (kind_cages)
        call seams_cages_destroy(box%raw)
      case (kind_topo)
        call seams_fingerprint_destroy(box%raw)
      case (kind_brine)
        call seams_brine_destroy(box%raw)
      end select
      box%raw = c_null_ptr
    end if
    deallocate(box)
    box => null()
  end subroutine handle_release

  !> Share `rhs` with `lhs`, or release `lhs` when `rhs` is disassociated.
  subroutine handle_share(lhs, rhs)
    type(handle_box), pointer, intent(inout) :: lhs
    type(handle_box), pointer, intent(in) :: rhs
    if (associated(lhs, rhs)) return
    call handle_release(lhs)
    if (associated(rhs)) then
      lhs => rhs
      lhs%refcount = lhs%refcount + 1
    else
      lhs => null()
    end if
  end subroutine handle_share

  !> The C pointer, or `c_null_ptr` when the box is closed.
  function handle_get(box) result(raw)
    type(handle_box), pointer, intent(in) :: box
    type(c_ptr) :: raw
    if (.not. associated(box)) then
      raw = c_null_ptr
    else
      raw = box%raw
    end if
  end function handle_get

  !> Copy a C string. A null pointer becomes an empty Fortran string.
  function c_to_f(ptr) result(text)
    type(c_ptr), intent(in) :: ptr
    character(len=:), allocatable :: text
    character(kind=c_char), pointer :: chars(:)
    integer :: i, n
    if (.not. c_associated(ptr)) then
      text = ""
      return
    end if
    n = int(strlen(ptr))
    allocate(character(len=n) :: text)
    if (n == 0) return
    call c_f_pointer(ptr, chars, [n])
    do i = 1, n
      text(i:i) = chars(i)
    end do
  end function c_to_f

  !> NUL-terminated copy of a Fortran string.
  function f_to_c(text) result(buf)
    character(len=*), intent(in) :: text
    character(kind=c_char), allocatable :: buf(:)
    integer :: i, n
    n = len_trim(text)
    allocate(buf(n + 1))
    do i = 1, n
      buf(i) = text(i:i)
    end do
    buf(n + 1) = c_null_char
  end function f_to_c

  !> Turn a queried C buffer into a Fortran string of length `n`.
  function copy_key(buf, n) result(text)
    character(kind=c_char), intent(in) :: buf(:)
    integer, intent(in) :: n
    character(len=:), allocatable :: text
    integer :: i
    allocate(character(len=n) :: text)
    do i = 1, n
      text(i:i) = buf(i)
    end do
  end function copy_key

  !> Intrinsic-style status. With neither `stat` nor `errmsg`, a failure stops.
  !!
  !! `code` 0 clears both arguments when they are present. A positive code
  !! stores the handle detail, or `detail` when that is present. The program
  !! stops only when the caller did not ask to receive the failure.
  subroutine check_status(code, handle, stat, errmsg, detail)
    integer(c_int), intent(in) :: code
    type(c_ptr), intent(in) :: handle
    integer(c_int), intent(out), optional :: stat
    character(len=:), allocatable, intent(out), optional :: errmsg
    character(len=*), intent(in), optional :: detail
    character(len=:), allocatable :: message
    if (code == SEAMS_OK) then
      if (present(stat)) stat = 0_c_int
      if (present(errmsg)) errmsg = ""
      return
    end if
    if (present(detail)) then
      message = detail
    else if (c_associated(handle)) then
      message = c_to_f(seams_handle_error(handle))
    else
      message = c_to_f(seams_last_error())
    end if
    if (len(message) == 0) message = "seams call failed"
    if (present(stat)) stat = code
    if (present(errmsg)) errmsg = message
    if (.not. present(stat) .and. .not. present(errmsg)) error stop message
  end subroutine check_status
end module seams_support
