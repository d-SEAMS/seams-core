! Copyright (c) d-SEAMS developers
!
!> Thin `iso_c_binding` image of `include/seams.h`.
!!
!! Every procedure is `bind(C)` to the C symbol of the same name,
!! except the string helpers at the bottom, which only pack Fortran
!! strings. Lengths that the C header types as `size_t` stay
!! `integer(c_size_t)`. A null pointer is `c_null_ptr`. This module
!! does not allocate result arrays and does not own handles: the
!! caller creates, sizes, and destroys them. The `seams` module is
!! the Fortran API.
module seams_c
  use, intrinsic :: iso_c_binding
  implicit none
  !> The call succeeded.
  integer(c_int), parameter :: SEAMS_OK = 0
  !> A required pointer is null.
  integer(c_int), parameter :: SEAMS_ERR_NULL = 1
  !> A number is outside the accepted range.
  integer(c_int), parameter :: SEAMS_ERR_RANGE = 2
  !> The cell is missing or has no volume.
  integer(c_int), parameter :: SEAMS_ERR_BOX = 3
  !> A caller argument is the wrong shape.
  integer(c_int), parameter :: SEAMS_ERR_ARG = 4
  !> The handle has no frame yet.
  integer(c_int), parameter :: SEAMS_ERR_STATE = 5
  !> An allocation failed.
  integer(c_int), parameter :: SEAMS_ERR_NOMEM = 6
  !> The library hit an internal failure.
  integer(c_int), parameter :: SEAMS_ERR_INTERNAL = 7
  !> CHILL+ rule from Nguyen and Molinero, 2015. Inclusive c_ij edges.
  integer(c_int), parameter :: SEAMS_CHILL_PUBLISHED = 0
  !> CHILL+ rule used by OVITO. The same thresholds, open edges.
  integer(c_int), parameter :: SEAMS_CHILL_OVITO = 1
  !> Neighbour graph of every pair inside the cutoff.
  integer(c_int), parameter :: SEAMS_GRAPH_CUTOFF = 0
  !> Mutual four-nearest neighbour graph.
  integer(c_int), parameter :: SEAMS_GRAPH_MUTUAL = 1
  !> Cage graph: mutual seeds and union completions.
  integer(c_int), parameter :: SEAMS_GRAPH_SEEDED = 2
  !> CHILL+ cubic-ice label.
  integer(c_int), parameter :: SEAMS_LABEL_CUBIC = 0

  interface
  !> Version text, for example 2.14.0. The pointer is static; the caller does not free it.
    function seams_version_string() bind(C, name="seams_version_string")
      import :: c_ptr
      type(c_ptr) :: seams_version_string
    end function

  !> Thread-local detail of the last failing call on this thread. Never NULL.
    function seams_last_error() bind(C, name="seams_last_error")
      import :: c_ptr
      type(c_ptr) :: seams_last_error
    end function

  !> Set the OpenMP team used by later analysis. Zero restores the saved default.
    function seams_set_num_threads(n) bind(C, name="seams_set_num_threads")
      import :: c_int
      integer(c_int), value :: n
      integer(c_int) :: seams_set_num_threads
    end function

  !> Team width the next analysis call will request.
    function seams_get_num_threads() bind(C, name="seams_get_num_threads")
      import :: c_int
      integer(c_int) :: seams_get_num_threads
    end function

  !> Allocate an empty frame. NULL means allocation failed.
    function seams_system_create() bind(C, name="seams_system_create")
      import :: c_ptr
      type(c_ptr) :: seams_system_create
    end function

  !> Free a frame. NULL is ignored. Neighbour and ring handles keep their own copies.
    subroutine seams_system_destroy(sys) bind(C, name="seams_system_destroy")
      import :: c_ptr
      type(c_ptr), value :: sys
    end subroutine

  !> Store points and an orthorhombic cell. types and ids may be c_null_ptr.
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

  !> Store points and a nine-number LAMMPS cell.
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

  !> CHILL+ labels and nine population counts. counts is label order, 1-based in Fortran.
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

  !> Local q4 and q6 on the cutoff shell.
    function seams_q4q6(sys, cutoff, typ, q4, q6) bind(C, name="seams_q4q6")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: sys
      real(c_double), value :: cutoff
      integer(c_int), value :: typ
      real(c_double), intent(out) :: q4(*)
      real(c_double), intent(out) :: q6(*)
      integer(c_int) :: seams_q4q6
    end function

  !> Build a cutoff or mutual neighbour list. out receives the handle.
    function seams_neighbours_build(sys, cutoff, typ, graph, out) &
        bind(C, name="seams_neighbours_build")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: sys
      real(c_double), value :: cutoff
      integer(c_int), value :: typ, graph
      type(c_ptr), intent(out) :: out
      integer(c_int) :: seams_neighbours_build
    end function

  !> Free a neighbour list. NULL is ignored.
    subroutine seams_neighbours_destroy(nb) bind(C, name="seams_neighbours_destroy")
      import :: c_ptr
      type(c_ptr), value :: nb
    end subroutine

  !> Borrow the neighbour CSR. The pointers are owned by the handle.
    function seams_neighbours_csr(nb, offsets, n_offsets, cols, nnz) &
        bind(C, name="seams_neighbours_csr")
      import :: c_ptr, c_int
      type(c_ptr), value :: nb
      type(c_ptr), intent(out) :: offsets, cols
      integer(c_int), intent(out) :: n_offsets, nnz
      integer(c_int) :: seams_neighbours_csr
    end function

  !> Primitive rings up to max_depth. out receives the handle.
    function seams_rings_find(sys, cutoff, typ, graph, max_depth, out) &
        bind(C, name="seams_rings_find")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: sys
      real(c_double), value :: cutoff
      integer(c_int), value :: typ, graph, max_depth
      type(c_ptr), intent(out) :: out
      integer(c_int) :: seams_rings_find
    end function

  !> Free a ring set. NULL is ignored.
    subroutine seams_rings_destroy(rings) bind(C, name="seams_rings_destroy")
      import :: c_ptr
      type(c_ptr), value :: rings
    end subroutine

  !> Count rings by size. counts(k) is the number of rings of size k-1 when the array is 1-based; a 0-based array stores size k at index k.
    function seams_rings_histogram(rings, counts, counts_len) &
        bind(C, name="seams_rings_histogram")
      import :: c_ptr, c_int
      type(c_ptr), value :: rings
      integer(c_int), intent(out) :: counts(*)
      integer(c_int), value :: counts_len
      integer(c_int) :: seams_rings_histogram
    end function

  !> Fraction of atoms that sit on a ring. size 0 counts every ring.
    function seams_rings_fraction(rings, size, fraction) bind(C, name="seams_rings_fraction")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: rings
      integer(c_int), value :: size
      real(c_double), intent(out) :: fraction
      integer(c_int) :: seams_rings_fraction
    end function

  !> Detail string stored on a handle. Survives a later call on another handle.
    function seams_handle_error(handle) bind(C, name="seams_handle_error")
      import :: c_ptr
      type(c_ptr), value :: handle
      type(c_ptr) :: seams_handle_error
    end function

  !> Number of points in the frame, or 0 when the handle is empty.
    function seams_system_count(sys) bind(C, name="seams_system_count")
      import :: c_ptr, c_int
      type(c_ptr), value :: sys
      integer(c_int) :: seams_system_count
    end function

  !> Read one 1-based frame. type_filter negative keeps every atom.
    function seams_system_read(sys, path, frame, type_filter) bind(C, name="seams_system_read")
      import :: c_ptr, c_int, c_char
      type(c_ptr), value :: sys
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int), value :: frame, type_filter
      integer(c_int) :: seams_system_read
    end function

  !> Original CHILL labels. labels_len is at least the frame count.
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

  !> Local q_l and the Lechner-Dellago average on the cutoff shell.
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

  !> Hexagonal and double-diamond cage membership.
    function seams_cages(sys, cutoff, typ, graph, k, complete, out) bind(C, name="seams_cages")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: sys
      real(c_double), value :: cutoff
      integer(c_int), value :: typ, graph, k, complete
      type(c_ptr), intent(out) :: out
      integer(c_int) :: seams_cages
    end function

  !> Free a cage-membership handle. NULL is ignored.
    subroutine seams_cages_destroy(cages) bind(C, name="seams_cages_destroy")
      import :: c_ptr
      type(c_ptr), value :: cages
    end subroutine

  !> Disjoint cage population counts. The four counts sum to the atoms the search kept.
    function seams_cages_counts(cages, n_hex, n_cubic, n_water, n_both) &
        bind(C, name="seams_cages_counts")
      import :: c_ptr, c_int
      type(c_ptr), value :: cages
      integer(c_int), intent(out) :: n_hex, n_cubic, n_water, n_both
      integer(c_int) :: seams_cages_counts
    end function

  !> Closed cages whose ring census matches spec, for example 512.
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

  !> k-nearest neighbour list. mutual nonzero keeps the mutual graph.
    function seams_neighbours_knn(sys, k, cutoff, typ, mutual, out) &
        bind(C, name="seams_neighbours_knn")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: sys
      integer(c_int), value :: k, typ, mutual
      real(c_double), value :: cutoff
      type(c_ptr), intent(out) :: out
      integer(c_int) :: seams_neighbours_knn
    end function

  !> Partial radial distribution. A null column queries the length.
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

  !> Running coordination number along the partial RDF.
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

  !> Coordination number integrated to rmax.
    function seams_coordination(sys, rmax, bins, type_i, type_j, cn) &
        bind(C, name="seams_coordination")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: sys
      real(c_double), value :: rmax
      integer(c_int), value :: bins, type_i, type_j
      real(c_double), intent(out) :: cn
      integer(c_int) :: seams_coordination
    end function

  !> Hydrogen-bond count. dist_cutoff is angstrom and angle_cutoff is degrees.
    function seams_hbonds(sys, cutoff, oxygen_type, hydrogen_type, dist_cutoff, angle_cutoff, &
        n_bonds) bind(C, name="seams_hbonds")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: sys
      real(c_double), value :: cutoff, dist_cutoff, angle_cutoff
      integer(c_int), value :: oxygen_type, hydrogen_type
      integer(c_int), intent(out) :: n_bonds
      integer(c_int) :: seams_hbonds
    end function

  !> Number density along one cell axis. axis is 0, 1, or 2.
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

  !> Topological fingerprint of the bonded graph.
    function seams_fingerprint(sys, cutoff, typ, k, graph, hops, colour_types, out) &
        bind(C, name="seams_fingerprint")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: sys
      real(c_double), value :: cutoff
      integer(c_int), value :: typ, k, graph, hops, colour_types
      type(c_ptr), intent(out) :: out
      integer(c_int) :: seams_fingerprint
    end function

  !> Free a fingerprint. NULL is ignored.
    subroutine seams_fingerprint_destroy(fp) bind(C, name="seams_fingerprint_destroy")
      import :: c_ptr
      type(c_ptr), value :: fp
    end subroutine

  !> Frame key text. A null buffer queries the length without the terminating NUL.
    function seams_fingerprint_key(fp, buf, cap, n_out) bind(C, name="seams_fingerprint_key")
      import :: c_ptr, c_size_t, c_int
      type(c_ptr), value :: fp, buf
      integer(c_size_t), value :: cap
      integer(c_size_t), intent(out) :: n_out
      integer(c_int) :: seams_fingerprint_key
    end function

  !> Primitive-ring census stored on the fingerprint. counts(k) is size k for a 0-based buffer.
    function seams_fingerprint_rings(fp, counts, counts_len) &
        bind(C, name="seams_fingerprint_rings")
      import :: c_ptr, c_int, c_size_t
      type(c_ptr), value :: fp
      integer(c_int), intent(out) :: counts(*)
      integer(c_size_t), value :: counts_len
      integer(c_int) :: seams_fingerprint_rings
    end function

  !> Mutual nearest unlike contacts of a site table.
    function seams_contact_pairs(sys, spec, n_pairs, n_cation, n_anion) &
        bind(C, name="seams_contact_pairs")
      import :: c_ptr, c_char, c_size_t, c_int
      type(c_ptr), value :: sys
      character(kind=c_char), intent(in) :: spec(*)
      integer(c_size_t), intent(out) :: n_pairs
      integer(c_int), intent(out) :: n_cation, n_anion
      integer(c_int) :: seams_contact_pairs
    end function

  !> Largest polar or apolar domain. subset 0 is apolar and 1 is polar.
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

  !> Perfect and imperfect prism counts by ring size.
    function seams_prisms(sys, cutoff, typ, perfect, imperfect, len) bind(C, name="seams_prisms")
      import :: c_ptr, c_double, c_int, c_size_t
      type(c_ptr), value :: sys
      real(c_double), value :: cutoff
      integer(c_int), value :: typ
      integer(c_int), intent(out) :: perfect(*), imperfect(*)
      integer(c_size_t), value :: len
      integer(c_int) :: seams_prisms
    end function

  !> Primitive-ring histogram of a monolayer.
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

  !> Number of primitive rings.
    function seams_rings_count(rings) bind(C, name="seams_rings_count")
      import :: c_ptr, c_int
      type(c_ptr), value :: rings
      integer(c_int) :: seams_rings_count
    end function

  !> Atom indices of one ring. index is 0-based. A null buffer queries the length.
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

  !> Per-atom site kind from a site-table spec.
    function seams_site_table(sys, spec, kinds, len) bind(C, name="seams_site_table")
      import :: c_ptr, c_char, c_int, c_size_t
      type(c_ptr), value :: sys
      character(kind=c_char), intent(in) :: spec(*)
      integer(c_int), intent(out) :: kinds(*)
      integer(c_size_t), value :: len
      integer(c_int) :: seams_site_table
    end function

    !> ABI number compiled into this library. It is 1.
    function seams_abi() bind(C, name="seams_abi")
      import :: c_int
      integer(c_int) :: seams_abi
    end function

    !> Packed release number, (major << 16) | (minor << 8) | patch.
    function seams_version() bind(C, name="seams_version")
      import :: c_int
      integer(c_int) :: seams_version
    end function

    !> Major version.
    function seams_version_major() bind(C, name="seams_version_major")
      import :: c_int
      integer(c_int) :: seams_version_major
    end function

    !> Minor version.
    function seams_version_minor() bind(C, name="seams_version_minor")
      import :: c_int
      integer(c_int) :: seams_version_minor
    end function

    !> Patch version.
    function seams_version_patch() bind(C, name="seams_version_patch")
      import :: c_int
      integer(c_int) :: seams_version_patch
    end function

    !> Short name of a status code. The pointer is static.
    function seams_status_string(status) bind(C, name="seams_status_string")
      import :: c_int, c_ptr
      integer(c_int), value :: status
      type(c_ptr) :: seams_status_string
    end function

    !> Name of a CHILL label. Unknown codes return "unknown".
    function seams_label_name(label) bind(C, name="seams_label_name")
      import :: c_int, c_ptr
      integer(c_int), value :: label
      type(c_ptr) :: seams_label_name
    end function

    !> Name of a ring-symmetry code.
    function seams_ring_symmetry_name(symmetry) bind(C, name="seams_ring_symmetry_name")
      import :: c_int, c_ptr
      integer(c_int), value :: symmetry
      type(c_ptr) :: seams_ring_symmetry_name
    end function

    !> Store points and a cell. n and box_len are size_t. box_len must be 9.
    function seams_system_set_atoms(sys, xyz, n, box9, box_len, types, ids) &
        bind(C, name="seams_system_set_atoms")
      import :: c_ptr, c_double, c_size_t, c_int
      type(c_ptr), value :: sys
      real(c_double), intent(in) :: xyz(*)
      integer(c_size_t), value :: n
      real(c_double), intent(in) :: box9(*)
      integer(c_size_t), value :: box_len
      type(c_ptr), value :: types, ids
      integer(c_int) :: seams_system_set_atoms
    end function

    !> Number of frames in a trajectory path.
    function seams_path_nframes(path, nframes) bind(C, name="seams_path_nframes")
      import :: c_char, c_size_t, c_int
      character(kind=c_char), intent(in) :: path(*)
      integer(c_size_t), intent(out) :: nframes
      integer(c_int) :: seams_path_nframes
    end function

    !> Number of neighbour rows.
    function seams_neighbours_count(nb) bind(C, name="seams_neighbours_count")
      import :: c_ptr, c_int
      type(c_ptr), value :: nb
      integer(c_int) :: seams_neighbours_count
    end function

    !> Copy the neighbour CSR into caller buffers. A null buffer queries the length.
    function seams_neighbours_copy(nb, offsets, offsets_cap, n_offsets, cols, cols_cap, nnz) &
        bind(C, name="seams_neighbours_copy")
      import :: c_ptr, c_size_t, c_int
      type(c_ptr), value :: nb, offsets, cols
      integer(c_size_t), value :: offsets_cap, cols_cap
      integer(c_size_t), intent(out) :: n_offsets, nnz
      integer(c_int) :: seams_neighbours_copy
    end function

    !> Local q_l and, on request, its neighbour average.
    function seams_ql(sys, cutoff, typ, order_l, graph, ql, ql_len, qlbar, qlbar_len) &
        bind(C, name="seams_ql")
      import :: c_ptr, c_double, c_int, c_size_t
      type(c_ptr), value :: sys, ql, qlbar
      real(c_double), value :: cutoff
      integer(c_int), value :: typ, order_l, graph
      integer(c_size_t), value :: ql_len, qlbar_len
      integer(c_int) :: seams_ql
    end function

    !> Lechner-Dellago q6-bar and the ice count above ice_cut.
    function seams_q6bar(sys, cutoff, typ, ice_cut, q6bar, q6bar_len, mean, median, n_ice) &
        bind(C, name="seams_q6bar")
      import :: c_ptr, c_double, c_int, c_size_t
      type(c_ptr), value :: sys, q6bar, mean, median, n_ice
      real(c_double), value :: cutoff, ice_cut
      integer(c_int), value :: typ
      integer(c_size_t), value :: q6bar_len
      integer(c_int) :: seams_q6bar
    end function

    !> Size and symmetry of each ring. Either pointer may be c_null_ptr, not both.
    function seams_rings_labels(rings, sizes, symmetry) bind(C, name="seams_rings_labels")
      import :: c_ptr, c_int
      type(c_ptr), value :: rings, sizes, symmetry
      integer(c_int) :: seams_rings_labels
    end function

    !> Per-atom cage code. len is at least the frame count stored on the handle.
    function seams_cages_labels(cages, per_atom, len) bind(C, name="seams_cages_labels")
      import :: c_ptr, c_int, c_size_t
      type(c_ptr), value :: cages
      integer(c_int), intent(out) :: per_atom(*)
      integer(c_size_t), value :: len
      integer(c_int) :: seams_cages_labels
    end function

    !> Per-atom ice shell label. A non-positive shell_rmsd selects the library ceiling.
    function seams_system_symmetry(sys, shell_rmsd, labels, labels_len) &
        bind(C, name="seams_system_symmetry")
      import :: c_ptr, c_double, c_int, c_size_t
      type(c_ptr), value :: sys
      real(c_double), value :: shell_rmsd
      integer(c_int), intent(out) :: labels(*)
      integer(c_size_t), value :: labels_len
      integer(c_int) :: seams_system_symmetry
    end function

    !> Planar ice front. ice(i) nonzero means atom i is ice.
    function seams_brine_analyse(sys, ice, n_ice, water_type, ion_type, n_bins, axis, out) &
        bind(C, name="seams_brine_analyse")
      import :: c_ptr, c_int, c_size_t
      type(c_ptr), value :: sys
      integer(c_int), intent(in) :: ice(*)
      integer(c_size_t), value :: n_ice
      integer(c_int), value :: water_type, ion_type, n_bins, axis
      type(c_ptr), intent(out) :: out
      integer(c_int) :: seams_brine_analyse
    end function

    !> Free a brine-front handle. NULL is ignored.
    subroutine seams_brine_destroy(front) bind(C, name="seams_brine_destroy")
      import :: c_ptr
      type(c_ptr), value :: front
    end subroutine

    !> The two planar crossings of the ice slab. Lengths are angstrom.
    function seams_brine_fronts(front, z_enter, z_exit, thickness, found) &
        bind(C, name="seams_brine_fronts")
      import :: c_ptr, c_int
      type(c_ptr), value :: front, z_enter, z_exit, thickness, found
      integer(c_int) :: seams_brine_fronts
    end function

    !> Molal rejection coefficient k = c_ice / c_liquid.
    function seams_brine_rejection(front, c_ice, c_liquid, k, defined) &
        bind(C, name="seams_brine_rejection")
      import :: c_ptr, c_int
      type(c_ptr), value :: front, c_ice, c_liquid, k, defined
      integer(c_int) :: seams_brine_rejection
    end function

    !> Concentration profile. A null column is skipped. cap is the capacity of each column.
    function seams_brine_profile(front, z_mid, molality, ice_fraction, cap, n_out) &
        bind(C, name="seams_brine_profile")
      import :: c_ptr, c_size_t, c_int
      type(c_ptr), value :: front, z_mid, molality, ice_fraction
      integer(c_size_t), value :: cap
      integer(c_size_t), intent(out) :: n_out
      integer(c_int) :: seams_brine_profile
    end function

    !> Brine pockets trapped in the ice, excluding the liquid reservoir.
    function seams_brine_pockets(front, n_pockets, reservoir_waters) &
        bind(C, name="seams_brine_pockets")
      import :: c_ptr, c_int
      type(c_ptr), value :: front, n_pockets, reservoir_waters
      integer(c_int) :: seams_brine_pockets
    end function

    !> Local key of one atom. Same buffer rule as seams_fingerprint_key.
    function seams_fingerprint_atom(fp, index, buf, cap, n_out) &
        bind(C, name="seams_fingerprint_atom")
      import :: c_ptr, c_int, c_size_t
      type(c_ptr), value :: fp, buf
      integer(c_int), value :: index
      integer(c_size_t), value :: cap
      integer(c_size_t), intent(out) :: n_out
      integer(c_int) :: seams_fingerprint_atom
    end function

  end interface

contains

  !> Copy a Fortran string into a NUL-terminated C buffer.
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

  !> Read one frame. path is a Fortran string.
  function read_frame(sys, path, frame, type_filter) result(code)
    type(c_ptr), value :: sys
    character(len=*), intent(in) :: path
    integer, intent(in) :: frame, type_filter
    integer(c_int) :: code
    character(kind=c_char), allocatable :: text(:)
    text = c_chars(path)
    code = seams_system_read(sys, text, int(frame, c_int), int(type_filter, c_int))
  end function

  !> Count cages of spec on the cutoff graph at 3.5 angstrom.
  function signature_count(sys, spec, n_cages, n_atoms) result(code)
    type(c_ptr), value :: sys
    character(len=*), intent(in) :: spec
    integer(c_size_t), intent(out) :: n_cages, n_atoms
    integer(c_int) :: code
    character(kind=c_char), allocatable :: text(:)
    text = c_chars(spec)
    code = seams_signature_cages(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, text, n_cages, n_atoms)
  end function

  !> Per-atom site kinds for spec.
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

  !> Contact-pair census for spec.
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

  !> Apolar domain size for spec at 3.5 angstrom.
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
end module seams_c
