! Copyright (c) d-SEAMS developers
module seams
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int), parameter :: SEAMS_OK = 0
  integer(c_int), parameter :: SEAMS_CHILL_PUBLISHED = 0
  integer(c_int), parameter :: SEAMS_CHILL_OVITO = 1
  integer(c_int), parameter :: SEAMS_GRAPH_CUTOFF = 0
  integer(c_int), parameter :: SEAMS_GRAPH_MUTUAL = 1
  integer(c_int), parameter :: SEAMS_GRAPH_SEEDED = 2
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

    function seams_handle_error(handle) bind(C, name="seams_handle_error")
      import :: c_ptr
      type(c_ptr), value :: handle
      type(c_ptr) :: seams_handle_error
    end function

    function seams_system_count(sys) bind(C, name="seams_system_count")
      import :: c_ptr, c_int
      type(c_ptr), value :: sys
      integer(c_int) :: seams_system_count
    end function

    function seams_system_read(sys, path, frame, type_filter) bind(C, name="seams_system_read")
      import :: c_ptr, c_int, c_char
      type(c_ptr), value :: sys
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int), value :: frame, type_filter
      integer(c_int) :: seams_system_read
    end function

    function seams_chill(sys, cutoff, typ, graph, rule, labels, labels_len, counts) &
        bind(C, name="seams_chill")
      import :: c_ptr, c_int, c_double, c_size_t
      type(c_ptr), value :: sys
      real(c_double), value :: cutoff
      integer(c_int), value :: typ, graph, rule
      integer(c_int), intent(out) :: labels(*)
      integer(c_size_t), value :: labels_len
      integer(c_int), intent(out) :: counts(*)
      integer(c_int) :: seams_chill
    end function

    function seams_steinhardt(sys, cutoff, typ, order_l, ql, ql_len, qlbar, qlbar_len) &
        bind(C, name="seams_steinhardt")
      import :: c_ptr, c_int, c_double, c_size_t
      type(c_ptr), value :: sys
      real(c_double), value :: cutoff
      integer(c_int), value :: typ, order_l
      real(c_double), intent(out) :: ql(*), qlbar(*)
      integer(c_size_t), value :: ql_len, qlbar_len
      integer(c_int) :: seams_steinhardt
    end function

    function seams_cages(sys, cutoff, typ, graph, k, complete, out) bind(C, name="seams_cages")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: sys
      real(c_double), value :: cutoff
      integer(c_int), value :: typ, graph, k, complete
      type(c_ptr), intent(out) :: out
      integer(c_int) :: seams_cages
    end function

    subroutine seams_cages_destroy(cages) bind(C, name="seams_cages_destroy")
      import :: c_ptr
      type(c_ptr), value :: cages
    end subroutine

    function seams_cages_counts(cages, n_hex, n_cubic, n_water, n_both) &
        bind(C, name="seams_cages_counts")
      import :: c_ptr, c_int
      type(c_ptr), value :: cages
      integer(c_int), intent(out) :: n_hex, n_cubic, n_water, n_both
      integer(c_int) :: seams_cages_counts
    end function

    function seams_signature_cages(sys, cutoff, typ, graph, spec, n_cages, n_atoms) &
        bind(C, name="seams_signature_cages")
      import :: c_ptr, c_int, c_double, c_char, c_size_t
      type(c_ptr), value :: sys
      real(c_double), value :: cutoff
      integer(c_int), value :: typ, graph
      character(kind=c_char), intent(in) :: spec(*)
      integer(c_size_t), intent(out) :: n_cages, n_atoms
      integer(c_int) :: seams_signature_cages
    end function

    function seams_neighbours_knn(sys, k, cutoff, typ, mutual, out) &
        bind(C, name="seams_neighbours_knn")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: sys
      integer(c_int), value :: k, typ, mutual
      real(c_double), value :: cutoff
      type(c_ptr), intent(out) :: out
      integer(c_int) :: seams_neighbours_knn
    end function

    function seams_rdf(sys, rmax, bins, type_i, type_j, r, g, count, cap, n_out, volume) &
        bind(C, name="seams_rdf")
      import :: c_ptr, c_int, c_double, c_size_t
      type(c_ptr), value :: sys
      real(c_double), value :: rmax
      integer(c_int), value :: bins, type_i, type_j
      type(c_ptr), value :: r, g, count, volume
      integer(c_size_t), value :: cap
      integer(c_size_t), intent(out) :: n_out
      integer(c_int) :: seams_rdf
    end function

    function seams_running_cn(sys, rmax, bins, type_i, type_j, r, cn, cap, n_out) &
        bind(C, name="seams_running_cn")
      import :: c_ptr, c_int, c_double, c_size_t
      type(c_ptr), value :: sys
      real(c_double), value :: rmax
      integer(c_int), value :: bins, type_i, type_j
      type(c_ptr), value :: r, cn
      integer(c_size_t), value :: cap
      integer(c_size_t), intent(out) :: n_out
      integer(c_int) :: seams_running_cn
    end function

    function seams_coordination(sys, rmax, bins, type_i, type_j, cn) &
        bind(C, name="seams_coordination")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: sys
      real(c_double), value :: rmax
      integer(c_int), value :: bins, type_i, type_j
      real(c_double), intent(out) :: cn
      integer(c_int) :: seams_coordination
    end function

    function seams_hbonds(sys, cutoff, oxygen_type, hydrogen_type, dist_cutoff, angle_cutoff, &
        n_bonds) bind(C, name="seams_hbonds")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: sys
      real(c_double), value :: cutoff, dist_cutoff, angle_cutoff
      integer(c_int), value :: oxygen_type, hydrogen_type
      integer(c_int), intent(out) :: n_bonds
      integer(c_int) :: seams_hbonds
    end function

    function seams_density_axis(sys, typ, bins, axis, z, rho, cap, n_out) &
        bind(C, name="seams_density_axis")
      import :: c_ptr, c_int, c_double, c_size_t
      type(c_ptr), value :: sys
      integer(c_int), value :: typ, bins, axis
      type(c_ptr), value :: z, rho
      integer(c_size_t), value :: cap
      integer(c_size_t), intent(out) :: n_out
      integer(c_int) :: seams_density_axis
    end function

    function seams_fingerprint(sys, cutoff, typ, k, graph, hops, colour_types, out) &
        bind(C, name="seams_fingerprint")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: sys
      real(c_double), value :: cutoff
      integer(c_int), value :: typ, k, graph, hops, colour_types
      type(c_ptr), intent(out) :: out
      integer(c_int) :: seams_fingerprint
    end function

    subroutine seams_fingerprint_destroy(fp) bind(C, name="seams_fingerprint_destroy")
      import :: c_ptr
      type(c_ptr), value :: fp
    end subroutine

    function seams_fingerprint_key(fp, buf, cap, n_out) bind(C, name="seams_fingerprint_key")
      import :: c_ptr, c_size_t, c_int
      type(c_ptr), value :: fp, buf
      integer(c_size_t), value :: cap
      integer(c_size_t), intent(out) :: n_out
      integer(c_int) :: seams_fingerprint_key
    end function

    function seams_fingerprint_rings(fp, counts, counts_len) &
        bind(C, name="seams_fingerprint_rings")
      import :: c_ptr, c_int, c_size_t
      type(c_ptr), value :: fp
      integer(c_int), intent(out) :: counts(*)
      integer(c_size_t), value :: counts_len
      integer(c_int) :: seams_fingerprint_rings
    end function

    function seams_contact_pairs(sys, spec, n_pairs, n_cation, n_anion) &
        bind(C, name="seams_contact_pairs")
      import :: c_ptr, c_char, c_size_t, c_int
      type(c_ptr), value :: sys
      character(kind=c_char), intent(in) :: spec(*)
      integer(c_size_t), intent(out) :: n_pairs
      integer(c_int), intent(out) :: n_cation, n_anion
      integer(c_int) :: seams_contact_pairs
    end function

    function seams_domains(sys, cutoff, spec, subset, n_subset, largest, percolation) &
        bind(C, name="seams_domains")
      import :: c_ptr, c_double, c_char, c_int
      type(c_ptr), value :: sys
      real(c_double), value :: cutoff
      character(kind=c_char), intent(in) :: spec(*)
      integer(c_int), value :: subset
      integer(c_int), intent(out) :: n_subset, largest
      real(c_double), intent(out) :: percolation
      integer(c_int) :: seams_domains
    end function

    function seams_prisms(sys, cutoff, typ, perfect, imperfect, len) bind(C, name="seams_prisms")
      import :: c_ptr, c_double, c_int, c_size_t
      type(c_ptr), value :: sys
      real(c_double), value :: cutoff
      integer(c_int), value :: typ
      integer(c_int), intent(out) :: perfect(*), imperfect(*)
      integer(c_size_t), value :: len
      integer(c_int) :: seams_prisms
    end function

    function seams_monolayer(sys, cutoff, typ, max_depth, sheet_area, counts, counts_len) &
        bind(C, name="seams_monolayer")
      import :: c_ptr, c_double, c_int, c_size_t
      type(c_ptr), value :: sys
      real(c_double), value :: cutoff, sheet_area
      integer(c_int), value :: typ, max_depth
      integer(c_int), intent(out) :: counts(*)
      integer(c_size_t), value :: counts_len
      integer(c_int) :: seams_monolayer
    end function

    function seams_rings_count(rings) bind(C, name="seams_rings_count")
      import :: c_ptr, c_int
      type(c_ptr), value :: rings
      integer(c_int) :: seams_rings_count
    end function

    function seams_rings_vertices(rings, index, atoms, cap, n_out) &
        bind(C, name="seams_rings_vertices")
      import :: c_ptr, c_int, c_size_t
      type(c_ptr), value :: rings
      integer(c_int), value :: index
      integer(c_int), intent(out) :: atoms(*)
      integer(c_size_t), value :: cap
      integer(c_size_t), intent(out) :: n_out
      integer(c_int) :: seams_rings_vertices
    end function

    function seams_site_table(sys, spec, kinds, len) bind(C, name="seams_site_table")
      import :: c_ptr, c_char, c_int, c_size_t
      type(c_ptr), value :: sys
      character(kind=c_char), intent(in) :: spec(*)
      integer(c_int), intent(out) :: kinds(*)
      integer(c_size_t), value :: len
      integer(c_int) :: seams_site_table
    end function
  end interface

contains

  function c_chars(text) result(buf)
    character(len=*), intent(in) :: text
    character(kind=c_char), allocatable :: buf(:)
    integer :: i, n
    n = len_trim(text)
    allocate(buf(n + 1))
    do i = 1, n
      buf(i) = text(i:i)
    end do
    buf(n + 1) = c_null_char
  end function

  function read_frame(sys, path, frame, type_filter) result(code)
    type(c_ptr), value :: sys
    character(len=*), intent(in) :: path
    integer, intent(in) :: frame, type_filter
    integer(c_int) :: code
    character(kind=c_char), allocatable :: text(:)
    text = c_chars(path)
    code = seams_system_read(sys, text, int(frame, c_int), int(type_filter, c_int))
  end function

  function signature_count(sys, spec, n_cages, n_atoms) result(code)
    type(c_ptr), value :: sys
    character(len=*), intent(in) :: spec
    integer(c_size_t), intent(out) :: n_cages, n_atoms
    integer(c_int) :: code
    character(kind=c_char), allocatable :: text(:)
    text = c_chars(spec)
    code = seams_signature_cages(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, text, n_cages, n_atoms)
  end function

  function site_kinds(sys, spec, kinds, n) result(code)
    type(c_ptr), value :: sys
    character(len=*), intent(in) :: spec
    integer(c_int), intent(out) :: kinds(*)
    integer(c_size_t), value :: n
    integer(c_int) :: code
    character(kind=c_char), allocatable :: text(:)
    text = c_chars(spec)
    code = seams_site_table(sys, text, kinds, n)
  end function

  function pair_count(sys, spec, n_pairs, n_cation, n_anion) result(code)
    type(c_ptr), value :: sys
    character(len=*), intent(in) :: spec
    integer(c_size_t), intent(out) :: n_pairs
    integer(c_int), intent(out) :: n_cation, n_anion
    integer(c_int) :: code
    character(kind=c_char), allocatable :: text(:)
    text = c_chars(spec)
    code = seams_contact_pairs(sys, text, n_pairs, n_cation, n_anion)
  end function

  function domain_size(sys, spec, n_subset, largest, percolation) result(code)
    type(c_ptr), value :: sys
    character(len=*), intent(in) :: spec
    integer(c_int), intent(out) :: n_subset, largest
    real(c_double), intent(out) :: percolation
    integer(c_int) :: code
    character(kind=c_char), allocatable :: text(:)
    text = c_chars(spec)
    code = seams_domains(sys, 3.5d0, text, 0, n_subset, largest, percolation)
  end function
end module seams
