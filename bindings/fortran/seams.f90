! Copyright (c) d-SEAMS developers
module seams
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int), parameter :: SEAMS_OK = 0
  integer(c_int), parameter :: SEAMS_CHILL_PUBLISHED = 0
  integer(c_int), parameter :: SEAMS_CHILL_OVITO = 1
  integer(c_int), parameter :: SEAMS_GRAPH_CUTOFF = 0
  integer(c_int), parameter :: SEAMS_GRAPH_MUTUAL = 1
  integer(c_int), parameter :: SEAMS_LABEL_CUBIC = 0

  interface
    function seams_version_string() bind(C, name="seams_version_string")
      import :: c_ptr
      type(c_ptr) :: seams_version_string
    end function

    function seams_last_error() bind(C, name="seams_last_error")
      import :: c_ptr
      type(c_ptr) :: seams_last_error
    end function

    function seams_set_num_threads(n) bind(C, name="seams_set_num_threads")
      import :: c_int
      integer(c_int), value :: n
      integer(c_int) :: seams_set_num_threads
    end function

    function seams_get_num_threads() bind(C, name="seams_get_num_threads")
      import :: c_int
      integer(c_int) :: seams_get_num_threads
    end function

    function seams_system_create() bind(C, name="seams_system_create")
      import :: c_ptr
      type(c_ptr) :: seams_system_create
    end function

    subroutine seams_system_destroy(sys) bind(C, name="seams_system_destroy")
      import :: c_ptr
      type(c_ptr), value :: sys
    end subroutine

    function seams_system_set_ortho(sys, xyz, n, ox, oy, oz, lx, ly, lz, types, ids) &
        bind(C, name="seams_system_set_ortho")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: sys
      real(c_double), intent(in) :: xyz(*)
      integer(c_int), value :: n
      real(c_double), value :: ox, oy, oz, lx, ly, lz
      type(c_ptr), value :: types, ids
      integer(c_int) :: seams_system_set_ortho
    end function

    function seams_system_set_frame(sys, xyz, n, box9, types, ids) &
        bind(C, name="seams_system_set_frame")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: sys
      real(c_double), intent(in) :: xyz(*)
      integer(c_int), value :: n
      real(c_double), intent(in) :: box9(*)
      type(c_ptr), value :: types, ids
      integer(c_int) :: seams_system_set_frame
    end function

    function seams_chill_plus(sys, cutoff, typ, graph, rule, labels, counts) &
        bind(C, name="seams_chill_plus")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: sys
      real(c_double), value :: cutoff
      integer(c_int), value :: typ, graph, rule
      integer(c_int), intent(out) :: labels(*)
      integer(c_int), intent(out) :: counts(*)
      integer(c_int) :: seams_chill_plus
    end function

    function seams_q4q6(sys, cutoff, typ, q4, q6) bind(C, name="seams_q4q6")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: sys
      real(c_double), value :: cutoff
      integer(c_int), value :: typ
      real(c_double), intent(out) :: q4(*)
      real(c_double), intent(out) :: q6(*)
      integer(c_int) :: seams_q4q6
    end function

    function seams_neighbours_build(sys, cutoff, typ, graph, out) &
        bind(C, name="seams_neighbours_build")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: sys
      real(c_double), value :: cutoff
      integer(c_int), value :: typ, graph
      type(c_ptr), intent(out) :: out
      integer(c_int) :: seams_neighbours_build
    end function

    subroutine seams_neighbours_destroy(nb) bind(C, name="seams_neighbours_destroy")
      import :: c_ptr
      type(c_ptr), value :: nb
    end subroutine

    function seams_neighbours_csr(nb, offsets, n_offsets, cols, nnz) &
        bind(C, name="seams_neighbours_csr")
      import :: c_ptr, c_int
      type(c_ptr), value :: nb
      type(c_ptr), intent(out) :: offsets, cols
      integer(c_int), intent(out) :: n_offsets, nnz
      integer(c_int) :: seams_neighbours_csr
    end function

    function seams_rings_find(sys, cutoff, typ, graph, max_depth, out) &
        bind(C, name="seams_rings_find")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: sys
      real(c_double), value :: cutoff
      integer(c_int), value :: typ, graph, max_depth
      type(c_ptr), intent(out) :: out
      integer(c_int) :: seams_rings_find
    end function

    subroutine seams_rings_destroy(rings) bind(C, name="seams_rings_destroy")
      import :: c_ptr
      type(c_ptr), value :: rings
    end subroutine

    function seams_rings_histogram(rings, counts, counts_len) &
        bind(C, name="seams_rings_histogram")
      import :: c_ptr, c_int
      type(c_ptr), value :: rings
      integer(c_int), intent(out) :: counts(*)
      integer(c_int), value :: counts_len
      integer(c_int) :: seams_rings_histogram
    end function

    function seams_rings_fraction(rings, size, fraction) bind(C, name="seams_rings_fraction")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: rings
      integer(c_int), value :: size
      real(c_double), intent(out) :: fraction
      integer(c_int) :: seams_rings_fraction
    end function
  end interface
end module seams
