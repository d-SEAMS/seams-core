! Copyright (c) d-SEAMS developers
submodule (seams) seams_methods
  use seams_kinds, only: c_ptr, c_null_ptr, c_associated, c_loc, c_int, c_double, c_size_t, c_char
  use seams_c, only: seams_version_string, seams_label_name, seams_status_string, &
    seams_set_num_threads, seams_get_num_threads, seams_system_create, &
    seams_system_set_ortho, seams_system_set_frame, seams_system_read, seams_system_count, &
    seams_chill, seams_chill_plus, seams_q4q6, seams_ql, seams_steinhardt, seams_q6bar, &
    seams_neighbours_build, seams_neighbours_knn, seams_neighbours_count, seams_neighbours_copy, &
    seams_rings_find, seams_rings_count, seams_rings_histogram, seams_rings_fraction, &
    seams_rings_labels, seams_rings_vertices, c_seams_cages => seams_cages, seams_cages_counts, seams_cages_labels, &
    seams_signature_cages, seams_rdf, seams_running_cn, seams_coordination, seams_hbonds, &
    seams_density_axis, seams_fingerprint, seams_fingerprint_key, seams_fingerprint_atom, &
    seams_fingerprint_rings, seams_contact_pairs, seams_domains, seams_prisms, seams_monolayer, &
    seams_site_table, seams_system_symmetry, seams_brine_analyse, seams_brine_fronts, &
    seams_brine_rejection, seams_brine_profile, seams_brine_pockets, seams_path_nframes
  use seams_support, only: handle_acquire, handle_release, handle_share, handle_get, &
    c_to_f, f_to_c, check_status, copy_key, kind_system, kind_neighbours, kind_rings, &
    kind_cages, kind_topo, kind_brine
  implicit none
contains

  function opt_real(value, fallback) result(out)
    real(seams_real), intent(in), optional :: value
    real(seams_real), intent(in) :: fallback
    real(seams_real) :: out
    if (present(value)) then
      out = value
    else
      out = fallback
    end if
  end function opt_real

  function opt_int(value, fallback) result(out)
    integer(seams_int), intent(in), optional :: value
    integer(seams_int), intent(in) :: fallback
    integer(c_int) :: out
    if (present(value)) then
      out = int(value, c_int)
    else
      out = int(fallback, c_int)
    end if
  end function opt_int

  function opt_flag(value, fallback) result(out)
    logical, intent(in), optional :: value
    logical, intent(in) :: fallback
    integer(c_int) :: out
    logical :: chosen
    if (present(value)) then
      chosen = value
    else
      chosen = fallback
    end if
    if (chosen) then
      out = 1_c_int
    else
      out = 0_c_int
    end if
  end function opt_flag

  function open_raw(box, stat, errmsg, what) result(raw)
    type(handle_box), pointer, intent(in) :: box
    integer(seams_int), intent(out), optional :: stat
    character(len=:), allocatable, intent(out), optional :: errmsg
    character(len=*), intent(in) :: what
    type(c_ptr) :: raw
    raw = handle_get(box)
    if (.not. c_associated(raw)) then
      call check_status(int(SEAMS_ERR_NULL, c_int), c_null_ptr, stat, errmsg, trim(what)//" is closed")
    end if
  end function open_raw

  function pack_optional(n, values, buf, ptr, stat, errmsg, what) result(ok)
    integer, intent(in) :: n
    integer(seams_int), intent(in), optional :: values(:)
    integer(c_int), allocatable, target, intent(out) :: buf(:)
    type(c_ptr), intent(out) :: ptr
    integer(seams_int), intent(out), optional :: stat
    character(len=:), allocatable, intent(out), optional :: errmsg
    character(len=*), intent(in) :: what
    logical :: ok
    ok = .true.
    ptr = c_null_ptr
    if (.not. present(values)) return
    if (size(values) /= n) then
      ok = .false.
      call check_status(int(SEAMS_ERR_ARG, c_int), c_null_ptr, stat, errmsg, trim(what)//" length is not the atom count")
      return
    end if
    allocate(buf(n))
    buf = int(values, c_int)
    ptr = c_loc(buf)
  end function pack_optional

  function shape_ok(xyz, stat, errmsg) result(ok)
    real(seams_real), intent(in) :: xyz(:, :)
    integer(seams_int), intent(out), optional :: stat
    character(len=:), allocatable, intent(out), optional :: errmsg
    logical :: ok
    ok = size(xyz, 1) == 3
    if (.not. ok) then
      call check_status(int(SEAMS_ERR_ARG, c_int), c_null_ptr, stat, errmsg, "xyz must have shape (3, n)")
    end if
  end function shape_ok

  subroutine take_key(raw, index, key, stat, errmsg)
    type(c_ptr), intent(in) :: raw
    integer(c_int), intent(in) :: index
    character(len=:), allocatable, intent(out) :: key
    integer(seams_int), intent(out), optional :: stat
    character(len=:), allocatable, intent(out), optional :: errmsg
    integer(c_size_t) :: n, cap
    integer(c_int) :: code
    character(kind=c_char), allocatable, target :: buf(:)
    if (index < 0_c_int) then
      code = seams_fingerprint_key(raw, c_null_ptr, 0_c_size_t, n)
    else
      code = seams_fingerprint_atom(raw, index, c_null_ptr, 0_c_size_t, n)
    end if
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      key = ""
      return
    end if
    allocate(buf(int(n) + 1))
    cap = n + 1_c_size_t
    if (index < 0_c_int) then
      code = seams_fingerprint_key(raw, c_loc(buf), cap, n)
    else
      code = seams_fingerprint_atom(raw, index, c_loc(buf), cap, n)
    end if
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      key = ""
      return
    end if
    key = copy_key(buf, int(n))
    call check_status(code, raw, stat, errmsg)
  end subroutine take_key

  module procedure new_frame
    type(c_ptr) :: raw
    raw = seams_system_create()
    if (.not. c_associated(raw)) then
      call check_status(int(SEAMS_ERR_NOMEM, c_int), c_null_ptr, stat, errmsg, "out of memory")
      return
    end if
    call handle_acquire(self%box, raw, kind_system)
    call check_status(0_c_int, raw, stat, errmsg)
  end procedure new_frame

  module procedure version_text
    text = c_to_f(seams_version_string())
  end procedure version_text

  module procedure set_num_threads
    integer(c_int) :: code
    code = seams_set_num_threads(int(n, c_int))
    call check_status(code, c_null_ptr, stat, errmsg)
  end procedure set_num_threads

  module procedure num_threads
    n = int(seams_get_num_threads(), seams_int)
  end procedure num_threads

  module procedure label_name
    text = c_to_f(seams_label_name(int(code, c_int)))
  end procedure label_name

  module procedure status_text
    text = c_to_f(seams_status_string(int(code, c_int)))
  end procedure status_text

  module procedure path_nframes
    character(kind=c_char), allocatable, target :: text(:)
    integer(c_size_t) :: count
    integer(c_int) :: code
    nframes = 0_seams_int
    text = f_to_c(path)
    code = seams_path_nframes(text, count)
    if (code /= SEAMS_OK) then
      call check_status(code, c_null_ptr, stat, errmsg)
      return
    end if
    nframes = int(count, seams_int)
    call check_status(code, c_null_ptr, stat, errmsg)
  end procedure path_nframes




  module procedure frame_set_ortho
    type(c_ptr) :: raw, types_ptr, ids_ptr
    integer(c_int), allocatable, target :: types_buf(:), ids_buf(:)
    integer(c_int) :: code, n
    real(seams_real), allocatable :: flat(:)
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    if (.not. shape_ok(xyz, stat, errmsg)) return
    n = int(size(xyz, 2), c_int)
    flat = reshape(xyz, [3 * size(xyz, 2)])
    if (.not. pack_optional(size(xyz, 2), types, types_buf, types_ptr, stat, errmsg, "types")) return
    if (.not. pack_optional(size(xyz, 2), ids, ids_buf, ids_ptr, stat, errmsg, "ids")) return
    code = seams_system_set_ortho(raw, flat, n, opt_real(ox, 0.0_seams_real), &
      opt_real(oy, 0.0_seams_real), opt_real(oz, 0.0_seams_real), lx, ly, lz, types_ptr, ids_ptr)
    call check_status(code, raw, stat, errmsg)
  end procedure frame_set_ortho

  module procedure frame_set_frame
    type(c_ptr) :: raw, types_ptr, ids_ptr
    integer(c_int), allocatable, target :: types_buf(:), ids_buf(:)
    integer(c_int) :: code, n
    real(seams_real), allocatable :: flat(:)
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    if (.not. shape_ok(xyz, stat, errmsg)) return
    n = int(size(xyz, 2), c_int)
    flat = reshape(xyz, [3 * size(xyz, 2)])
    if (.not. pack_optional(size(xyz, 2), types, types_buf, types_ptr, stat, errmsg, "types")) return
    if (.not. pack_optional(size(xyz, 2), ids, ids_buf, ids_ptr, stat, errmsg, "ids")) return
    code = seams_system_set_frame(raw, flat, n, box, types_ptr, ids_ptr)
    call check_status(code, raw, stat, errmsg)
  end procedure frame_set_frame

  module procedure frame_read
    type(c_ptr) :: raw
    character(kind=c_char), allocatable, target :: text(:)
    integer(c_int) :: code
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    text = f_to_c(path)
    code = seams_system_read(raw, text, opt_int(frame, 1_seams_int), opt_int(type_filter, -1_seams_int))
    call check_status(code, raw, stat, errmsg)
  end procedure frame_read

  module procedure frame_count
    type(c_ptr) :: raw
    raw = handle_get(self%box)
    if (.not. c_associated(raw)) then
      n = 0_seams_int
    else
      n = int(seams_system_count(raw), seams_int)
    end if
  end procedure frame_count

  subroutine chill_common(self, labels, counts, cutoff, typ, graph, rule, which, stat, errmsg)
    class(seams_frame), intent(inout) :: self
    integer(seams_int), allocatable, intent(out) :: labels(:)
    integer(seams_int), allocatable, intent(out), optional :: counts(:)
    real(seams_real), intent(in), optional :: cutoff
    integer(seams_int), intent(in), optional :: typ, graph, rule
    integer, intent(in) :: which
    integer(seams_int), intent(out), optional :: stat
    character(len=:), allocatable, intent(out), optional :: errmsg
    type(c_ptr) :: raw
    integer(c_int) :: code, n
    integer(seams_int), allocatable, target :: tally(:)
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    n = seams_system_count(raw)
    if (n <= 0_c_int) then
      call check_status(int(SEAMS_ERR_STATE, c_int), raw, stat, errmsg, "frame has no atoms")
      return
    end if
    allocate(labels(n), tally(9))
    if (which == 0) then
      code = seams_chill(raw, opt_real(cutoff, 3.5_seams_real), opt_int(typ, 0_seams_int), &
        opt_int(graph, SEAMS_GRAPH_CUTOFF), opt_int(rule, SEAMS_CHILL_PUBLISHED), labels, &
        int(n, c_size_t), tally)
    else
      code = seams_chill_plus(raw, opt_real(cutoff, 3.5_seams_real), opt_int(typ, 0_seams_int), &
        opt_int(graph, SEAMS_GRAPH_CUTOFF), opt_int(rule, SEAMS_CHILL_PUBLISHED), labels, tally)
    end if
    if (code /= SEAMS_OK) then
      deallocate(labels)
      call check_status(code, raw, stat, errmsg)
      return
    end if
    if (present(counts)) counts = tally
    call check_status(code, raw, stat, errmsg)
  end subroutine chill_common

  module procedure frame_chill
    call chill_common(self, labels, counts, cutoff, typ, graph, rule, 0, stat, errmsg)
  end procedure frame_chill

  module procedure frame_chill_plus
    call chill_common(self, labels, counts, cutoff, typ, graph, rule, 1, stat, errmsg)
  end procedure frame_chill_plus

  module procedure frame_q4q6
    type(c_ptr) :: raw
    integer(c_int) :: code, n
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    n = seams_system_count(raw)
    allocate(q4(n), q6(n))
    code = seams_q4q6(raw, opt_real(cutoff, 3.5_seams_real), opt_int(typ, 0_seams_int), q4, q6)
    if (code /= SEAMS_OK) then
      deallocate(q4, q6)
      call check_status(code, raw, stat, errmsg)
      return
    end if
    call check_status(code, raw, stat, errmsg)
  end procedure frame_q4q6

  module procedure frame_ql
    type(c_ptr) :: raw, bar_ptr
    integer(c_int) :: code, n
    integer(c_size_t) :: bar_len
    real(seams_real), allocatable, target :: qbuf(:), barbuf(:)
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    n = seams_system_count(raw)
    allocate(qbuf(n))
    bar_ptr = c_null_ptr
    bar_len = 0_c_size_t
    if (present(qlbar)) then
      allocate(barbuf(n))
      bar_ptr = c_loc(barbuf)
      bar_len = int(n, c_size_t)
    end if
    code = seams_ql(raw, opt_real(cutoff, 3.5_seams_real), opt_int(typ, 0_seams_int), &
      int(order_l, c_int), opt_int(graph, SEAMS_GRAPH_CUTOFF), c_loc(qbuf), int(n, c_size_t), &
      bar_ptr, bar_len)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    ql = qbuf
    if (present(qlbar)) qlbar = barbuf
    call check_status(code, raw, stat, errmsg)
  end procedure frame_ql

  module procedure frame_steinhardt
    type(c_ptr) :: raw
    integer(c_int) :: code, n
    real(seams_real), allocatable :: bar(:)
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    n = seams_system_count(raw)
    allocate(ql(n), bar(n))
    code = seams_steinhardt(raw, opt_real(cutoff, 3.5_seams_real), opt_int(typ, 0_seams_int), &
      int(order_l, c_int), ql, int(n, c_size_t), bar, int(n, c_size_t))
    if (code /= SEAMS_OK) then
      deallocate(ql)
      call check_status(code, raw, stat, errmsg)
      return
    end if
    if (present(qlbar)) qlbar = bar
    call check_status(code, raw, stat, errmsg)
  end procedure frame_steinhardt

  module procedure frame_q6bar
    type(c_ptr) :: raw, qptr, mean_ptr, median_ptr, ice_ptr
    integer(c_int) :: code, n
    integer(c_size_t) :: qlen
    real(seams_real), allocatable, target :: qbuf(:)
    real(seams_real), target :: mean_b, median_b
    integer(c_int), target :: ice_b
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    if (.not. (present(q6bar) .or. present(mean) .or. present(median) .or. present(n_ice))) then
      call check_status(int(SEAMS_ERR_ARG, c_int), raw, stat, errmsg, "q6bar needs one output")
      return
    end if
    n = seams_system_count(raw)
    qptr = c_null_ptr
    mean_ptr = c_null_ptr
    median_ptr = c_null_ptr
    ice_ptr = c_null_ptr
    qlen = 0_c_size_t
    if (present(q6bar)) then
      allocate(qbuf(n))
      qptr = c_loc(qbuf)
      qlen = int(n, c_size_t)
    end if
    if (present(mean)) mean_ptr = c_loc(mean_b)
    if (present(median)) median_ptr = c_loc(median_b)
    if (present(n_ice)) ice_ptr = c_loc(ice_b)
    code = seams_q6bar(raw, opt_real(cutoff, 3.5_seams_real), opt_int(typ, 0_seams_int), &
      opt_real(ice_cut, 0.36_seams_real), qptr, qlen, mean_ptr, median_ptr, ice_ptr)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    if (present(q6bar)) call move_alloc(qbuf, q6bar)
    if (present(mean)) mean = mean_b
    if (present(median)) median = median_b
    if (present(n_ice)) n_ice = int(ice_b, seams_int)
    call check_status(code, raw, stat, errmsg)
  end procedure frame_q6bar

  module procedure frame_neighbours
    type(c_ptr) :: raw, out
    integer(c_int) :: code
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    code = seams_neighbours_build(raw, opt_real(cutoff, 3.5_seams_real), opt_int(typ, 0_seams_int), &
      opt_int(graph, SEAMS_GRAPH_CUTOFF), out)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    call handle_acquire(list%box, out, kind_neighbours)
    call check_status(code, raw, stat, errmsg)
  end procedure frame_neighbours

  module procedure frame_knn
    type(c_ptr) :: raw, out
    integer(c_int) :: code
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    code = seams_neighbours_knn(raw, int(k, c_int), opt_real(cutoff, 3.5_seams_real), &
      opt_int(typ, 0_seams_int), opt_flag(mutual, .true.), out)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    call handle_acquire(list%box, out, kind_neighbours)
    call check_status(code, raw, stat, errmsg)
  end procedure frame_knn

  module procedure frame_rings
    type(c_ptr) :: raw, out
    integer(c_int) :: code, depth
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    depth = opt_int(max_depth, 6_seams_int)
    code = seams_rings_find(raw, opt_real(cutoff, 3.5_seams_real), opt_int(typ, 0_seams_int), &
      opt_int(graph, SEAMS_GRAPH_CUTOFF), depth, out)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    call handle_acquire(rings%box, out, kind_rings)
    rings%max_depth = int(depth, seams_int)
    call check_status(code, raw, stat, errmsg)
  end procedure frame_rings

  module procedure frame_cages
    type(c_ptr) :: raw, out
    integer(c_int) :: code
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    code = c_seams_cages(raw, opt_real(cutoff, 3.5_seams_real), opt_int(typ, 0_seams_int), &
      opt_int(graph, SEAMS_GRAPH_SEEDED), opt_int(k, 4_seams_int), opt_flag(complete, .false.), out)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    call handle_acquire(cages%box, out, kind_cages)
    cages%natoms = int(seams_system_count(raw), seams_int)
    call check_status(code, raw, stat, errmsg)
  end procedure frame_cages

  module procedure frame_signature
    type(c_ptr) :: raw
    character(kind=c_char), allocatable, target :: text(:)
    integer(c_size_t) :: cages_n, atoms_n
    integer(c_int) :: code
    n_cages = 0_seams_int
    n_atoms = 0_seams_int
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    text = f_to_c(spec)
    code = seams_signature_cages(raw, opt_real(cutoff, 3.5_seams_real), opt_int(typ, 0_seams_int), &
      opt_int(graph, SEAMS_GRAPH_CUTOFF), text, cages_n, atoms_n)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    n_cages = int(cages_n, seams_int)
    n_atoms = int(atoms_n, seams_int)
    call check_status(code, raw, stat, errmsg)
  end procedure frame_signature

  module procedure frame_rdf
    type(c_ptr) :: raw, count_ptr, volume_ptr
    integer(c_size_t) :: n, cap
    integer(c_int) :: code
    real(seams_real), allocatable, target :: rbuf(:), gbuf(:)
    integer(seams_int), allocatable, target :: cbuf(:)
    real(seams_real), target :: volume_b
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    code = seams_rdf(raw, opt_real(rmax, 6.0_seams_real), opt_int(bins, 0_seams_int), &
      opt_int(type_i, 0_seams_int), opt_int(type_j, 0_seams_int), c_null_ptr, c_null_ptr, &
      c_null_ptr, 0_c_size_t, n, c_null_ptr)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    allocate(rbuf(n), gbuf(n))
    count_ptr = c_null_ptr
    if (present(count)) then
      allocate(cbuf(n))
      count_ptr = c_loc(cbuf)
    end if
    volume_ptr = c_null_ptr
    if (present(volume)) volume_ptr = c_loc(volume_b)
    cap = n
    code = seams_rdf(raw, opt_real(rmax, 6.0_seams_real), opt_int(bins, 0_seams_int), &
      opt_int(type_i, 0_seams_int), opt_int(type_j, 0_seams_int), c_loc(rbuf), c_loc(gbuf), &
      count_ptr, cap, n, volume_ptr)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    r = rbuf
    g = gbuf
    if (present(count)) count = cbuf
    if (present(volume)) volume = volume_b
    call check_status(code, raw, stat, errmsg)
  end procedure frame_rdf

  module procedure frame_running_cn
    type(c_ptr) :: raw
    integer(c_size_t) :: n, cap
    integer(c_int) :: code
    real(seams_real), allocatable, target :: rbuf(:), cbuf(:)
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    code = seams_running_cn(raw, opt_real(rmax, 6.0_seams_real), opt_int(bins, 0_seams_int), &
      opt_int(type_i, 0_seams_int), opt_int(type_j, 0_seams_int), c_null_ptr, c_null_ptr, &
      0_c_size_t, n)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    allocate(rbuf(n), cbuf(n))
    cap = n
    code = seams_running_cn(raw, opt_real(rmax, 6.0_seams_real), opt_int(bins, 0_seams_int), &
      opt_int(type_i, 0_seams_int), opt_int(type_j, 0_seams_int), c_loc(rbuf), c_loc(cbuf), cap, n)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    r = rbuf
    cn = cbuf
    call check_status(code, raw, stat, errmsg)
  end procedure frame_running_cn

  module procedure frame_coordination
    type(c_ptr) :: raw
    integer(c_int) :: code
    real(c_double), target :: value
    cn = 0.0_seams_real
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    code = seams_coordination(raw, opt_real(rmax, 3.5_seams_real), opt_int(bins, 0_seams_int), &
      opt_int(type_i, 0_seams_int), opt_int(type_j, 0_seams_int), value)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    cn = value
    call check_status(code, raw, stat, errmsg)
  end procedure frame_coordination

  module procedure frame_hbonds
    type(c_ptr) :: raw
    integer(c_int) :: code, bonds
    n_bonds = 0_seams_int
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    code = seams_hbonds(raw, opt_real(cutoff, 3.5_seams_real), opt_int(oxygen_type, 1_seams_int), &
      opt_int(hydrogen_type, 2_seams_int), opt_real(dist, 2.42_seams_real), &
      opt_real(angle, 30.0_seams_real), bonds)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    n_bonds = int(bonds, seams_int)
    call check_status(code, raw, stat, errmsg)
  end procedure frame_hbonds

  module procedure frame_density
    type(c_ptr) :: raw
    integer(c_size_t) :: n, cap
    integer(c_int) :: code
    real(seams_real), allocatable, target :: zbuf(:), rbuf(:)
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    code = seams_density_axis(raw, opt_int(typ, 0_seams_int), opt_int(bins, 0_seams_int), &
      opt_int(axis, 2_seams_int), c_null_ptr, c_null_ptr, 0_c_size_t, n)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    allocate(zbuf(n), rbuf(n))
    cap = n
    code = seams_density_axis(raw, opt_int(typ, 0_seams_int), opt_int(bins, 0_seams_int), &
      opt_int(axis, 2_seams_int), c_loc(zbuf), c_loc(rbuf), cap, n)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    z = zbuf
    rho = rbuf
    call check_status(code, raw, stat, errmsg)
  end procedure frame_density

  module procedure frame_fingerprint
    type(c_ptr) :: raw, out
    integer(c_int) :: code
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    code = seams_fingerprint(raw, opt_real(cutoff, 3.5_seams_real), opt_int(typ, 0_seams_int), &
      opt_int(k, 4_seams_int), opt_int(graph, SEAMS_GRAPH_CUTOFF), opt_int(hops, 2_seams_int), &
      opt_flag(colours, .false.), out)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    call handle_acquire(fp%box, out, kind_topo)
    call check_status(code, raw, stat, errmsg)
  end procedure frame_fingerprint

  module procedure frame_contact_pairs
    type(c_ptr) :: raw
    character(kind=c_char), allocatable, target :: text(:)
    integer(c_size_t) :: pairs
    integer(c_int) :: code, cations, anions
    n_pairs = 0_seams_int
    n_cation = 0_seams_int
    n_anion = 0_seams_int
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    text = f_to_c(spec)
    code = seams_contact_pairs(raw, text, pairs, cations, anions)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    n_pairs = int(pairs, seams_int)
    n_cation = int(cations, seams_int)
    n_anion = int(anions, seams_int)
    call check_status(code, raw, stat, errmsg)
  end procedure frame_contact_pairs

  module procedure frame_domains
    type(c_ptr) :: raw
    character(kind=c_char), allocatable, target :: text(:)
    integer(c_int) :: code, subset_n, largest_n
    real(c_double), target :: perc
    n_subset = 0_seams_int
    largest = 0_seams_int
    percolation = 0.0_seams_real
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    text = f_to_c(spec)
    code = seams_domains(raw, opt_real(cutoff, 3.5_seams_real), text, opt_flag(polar, .false.), &
      subset_n, largest_n, perc)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    n_subset = int(subset_n, seams_int)
    largest = int(largest_n, seams_int)
    percolation = perc
    call check_status(code, raw, stat, errmsg)
  end procedure frame_domains

  module procedure frame_prisms
    type(c_ptr) :: raw
    integer(c_int) :: code
    integer(seams_int), allocatable, target :: perfect_b(:), imperfect_b(:)
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    allocate(perfect_b(0:6), imperfect_b(0:6))
    code = seams_prisms(raw, opt_real(cutoff, 3.5_seams_real), opt_int(typ, 0_seams_int), &
      perfect_b, imperfect_b, 7_c_size_t)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    perfect = perfect_b
    imperfect = imperfect_b
    call check_status(code, raw, stat, errmsg)
  end procedure frame_prisms

  module procedure frame_monolayer
    type(c_ptr) :: raw
    integer(c_int) :: code, depth
    integer(seams_int), allocatable, target :: buf(:)
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    depth = opt_int(max_depth, 6_seams_int)
    allocate(buf(0:depth))
    buf = 0_seams_int
    code = seams_monolayer(raw, opt_real(cutoff, 3.5_seams_real), opt_int(typ, 0_seams_int), &
      depth, sheet_area, buf, int(size(buf), c_size_t))
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    counts = buf
    call check_status(code, raw, stat, errmsg)
  end procedure frame_monolayer

  module procedure frame_site_kinds
    type(c_ptr) :: raw
    character(kind=c_char), allocatable, target :: text(:)
    integer(c_int) :: code, n
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    n = seams_system_count(raw)
    allocate(kinds(n))
    text = f_to_c(spec)
    code = seams_site_table(raw, text, kinds, int(n, c_size_t))
    if (code /= SEAMS_OK) then
      deallocate(kinds)
      call check_status(code, raw, stat, errmsg)
      return
    end if
    call check_status(code, raw, stat, errmsg)
  end procedure frame_site_kinds

  module procedure frame_symmetry
    type(c_ptr) :: raw
    integer(c_int) :: code, n
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    n = seams_system_count(raw)
    allocate(labels(n))
    code = seams_system_symmetry(raw, opt_real(shell_rmsd, 0.0_seams_real), labels, int(n, c_size_t))
    if (code /= SEAMS_OK) then
      deallocate(labels)
      call check_status(code, raw, stat, errmsg)
      return
    end if
    call check_status(code, raw, stat, errmsg)
  end procedure frame_symmetry

  module procedure frame_brine
    type(c_ptr) :: raw, out
    integer(c_int) :: code
    integer(c_int), allocatable, target :: flags(:)
    raw = open_raw(self%box, stat, errmsg, "frame")
    if (.not. c_associated(raw)) return
    if (size(ice) /= seams_system_count(raw)) then
      call check_status(int(SEAMS_ERR_ARG, c_int), raw, stat, errmsg, "ice length is not the atom count")
      return
    end if
    flags = int(ice, c_int)
    code = seams_brine_analyse(raw, flags, int(size(flags), c_size_t), opt_int(water_type, 0_seams_int), &
      opt_int(ion_type, 0_seams_int), opt_int(n_bins, 0_seams_int), opt_int(axis, 2_seams_int), out)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    call handle_acquire(front%box, out, kind_brine)
    call check_status(code, raw, stat, errmsg)
  end procedure frame_brine

  module procedure neighbour_count
    type(c_ptr) :: raw
    raw = handle_get(self%box)
    if (.not. c_associated(raw)) then
      n = 0_seams_int
    else
      n = int(seams_neighbours_count(raw), seams_int)
    end if
  end procedure neighbour_count

  subroutine neighbour_column(self, which, values, stat, errmsg)
    class(seams_neighbours), intent(in) :: self
    integer, intent(in) :: which
    integer(seams_int), allocatable, intent(out) :: values(:)
    integer(seams_int), intent(out), optional :: stat
    character(len=:), allocatable, intent(out), optional :: errmsg
    type(c_ptr) :: raw, off_ptr, col_ptr
    integer(c_size_t) :: n_off, nnz, off_cap, col_cap
    integer(c_int) :: code
    integer(c_int), allocatable, target :: off_buf(:), col_buf(:)
    raw = open_raw(self%box, stat, errmsg, "neighbours")
    if (.not. c_associated(raw)) return
    code = seams_neighbours_copy(raw, c_null_ptr, 0_c_size_t, n_off, c_null_ptr, 0_c_size_t, nnz)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    off_ptr = c_null_ptr
    col_ptr = c_null_ptr
    if (which == 0) then
      allocate(off_buf(n_off))
      off_ptr = c_loc(off_buf)
    else
      allocate(col_buf(nnz))
      col_ptr = c_loc(col_buf)
    end if
    off_cap = n_off
    col_cap = nnz
    code = seams_neighbours_copy(raw, off_ptr, off_cap, n_off, col_ptr, col_cap, nnz)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    if (which == 0) then
      values = int(off_buf, seams_int)
    else
      values = int(col_buf, seams_int)
    end if
    call check_status(code, raw, stat, errmsg)
  end subroutine neighbour_column

  module procedure neighbour_offsets
    call neighbour_column(self, 0, offsets, stat, errmsg)
  end procedure neighbour_offsets
  module procedure neighbour_columns
    call neighbour_column(self, 1, columns, stat, errmsg)
  end procedure neighbour_columns

  module procedure ring_count
    type(c_ptr) :: raw
    raw = handle_get(self%box)
    if (.not. c_associated(raw)) then
      n = 0_seams_int
    else
      n = int(seams_rings_count(raw), seams_int)
    end if
  end procedure ring_count
  module procedure ring_histogram
    type(c_ptr) :: raw
    integer(c_int) :: code
    integer(seams_int), allocatable, target :: buf(:)
    raw = open_raw(self%box, stat, errmsg, "rings")
    if (.not. c_associated(raw)) return
    allocate(buf(0:self%max_depth))
    buf = 0_seams_int
    code = seams_rings_histogram(raw, buf, int(size(buf), c_int))
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    counts = buf
    call check_status(code, raw, stat, errmsg)
  end procedure ring_histogram
  module procedure ring_fraction
    type(c_ptr) :: raw
    integer(c_int) :: code
    real(c_double), target :: value
    fraction = 0.0_seams_real
    raw = open_raw(self%box, stat, errmsg, "rings")
    if (.not. c_associated(raw)) return
    code = seams_rings_fraction(raw, int(size, c_int), value)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    fraction = value
    call check_status(code, raw, stat, errmsg)
  end procedure ring_fraction
  subroutine ring_one(self, which, values, stat, errmsg)
    class(seams_rings), intent(in) :: self
    integer, intent(in) :: which
    integer(seams_int), allocatable, intent(out) :: values(:)
    integer(seams_int), intent(out), optional :: stat
    character(len=:), allocatable, intent(out), optional :: errmsg
    type(c_ptr) :: raw, size_ptr, sym_ptr
    integer(c_int) :: code, n
    integer(c_int), allocatable, target :: size_buf(:), sym_buf(:)
    raw = open_raw(self%box, stat, errmsg, "rings")
    if (.not. c_associated(raw)) return
    n = seams_rings_count(raw)
    size_ptr = c_null_ptr
    sym_ptr = c_null_ptr
    if (which == 0) then
      allocate(size_buf(n))
      size_ptr = c_loc(size_buf)
    else
      allocate(sym_buf(n))
      sym_ptr = c_loc(sym_buf)
    end if
    code = seams_rings_labels(raw, size_ptr, sym_ptr)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    if (which == 0) then
      values = int(size_buf, seams_int)
    else
      values = int(sym_buf, seams_int)
    end if
    call check_status(code, raw, stat, errmsg)
  end subroutine ring_one
  module procedure ring_sizes
    call ring_one(self, 0, sizes, stat, errmsg)
  end procedure ring_sizes
  module procedure ring_symmetry
    call ring_one(self, 1, codes, stat, errmsg)
  end procedure ring_symmetry
  module procedure ring_vertices
    type(c_ptr) :: raw
    integer(c_int) :: code
    integer(c_size_t) :: n
    integer(c_int), allocatable, target :: buf(:)
    raw = open_raw(self%box, stat, errmsg, "rings")
    if (.not. c_associated(raw)) return
    allocate(buf(12))
    code = seams_rings_vertices(raw, int(index - 1, c_int), buf, 12_c_size_t, n)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    atoms = int(buf(1:int(n)), seams_int)
    call check_status(code, raw, stat, errmsg)
  end procedure ring_vertices

  module procedure cage_counts
    type(c_ptr) :: raw
    integer(c_int) :: code, hex, cubic, water, both
    n_hex = 0_seams_int
    n_cubic = 0_seams_int
    n_water = 0_seams_int
    n_both = 0_seams_int
    raw = open_raw(self%box, stat, errmsg, "cages")
    if (.not. c_associated(raw)) return
    code = seams_cages_counts(raw, hex, cubic, water, both)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    n_hex = int(hex, seams_int)
    n_cubic = int(cubic, seams_int)
    n_water = int(water, seams_int)
    n_both = int(both, seams_int)
    call check_status(code, raw, stat, errmsg)
  end procedure cage_counts
  module procedure cage_labels
    type(c_ptr) :: raw
    integer(c_int) :: code
    raw = open_raw(self%box, stat, errmsg, "cages")
    if (.not. c_associated(raw)) return
    allocate(labels(self%natoms))
    code = seams_cages_labels(raw, labels, int(self%natoms, c_size_t))
    if (code /= SEAMS_OK) then
      deallocate(labels)
      call check_status(code, raw, stat, errmsg)
      return
    end if
    call check_status(code, raw, stat, errmsg)
  end procedure cage_labels

  module procedure topo_key
    type(c_ptr) :: raw
    raw = open_raw(self%box, stat, errmsg, "fingerprint")
    if (.not. c_associated(raw)) return
    call take_key(raw, -1_c_int, key, stat, errmsg)
  end procedure topo_key
  module procedure topo_atom_key
    type(c_ptr) :: raw
    raw = open_raw(self%box, stat, errmsg, "fingerprint")
    if (.not. c_associated(raw)) return
    call take_key(raw, int(index - 1, c_int), key, stat, errmsg)
  end procedure topo_atom_key
  module procedure topo_rings
    type(c_ptr) :: raw
    integer(c_int) :: code
    integer(seams_int), allocatable, target :: buf(:)
    raw = open_raw(self%box, stat, errmsg, "fingerprint")
    if (.not. c_associated(raw)) return
    allocate(buf(0:12))
    buf = 0_seams_int
    code = seams_fingerprint_rings(raw, buf, int(size(buf), c_size_t))
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    counts = buf
    call check_status(code, raw, stat, errmsg)
  end procedure topo_rings

  module procedure brine_fronts
    type(c_ptr) :: raw, enter_ptr, exit_ptr, thick_ptr, found_ptr
    integer(c_int) :: code
    real(c_double), target :: enter_b, exit_b, thick_b
    integer(c_int), target :: found_b
    raw = open_raw(self%box, stat, errmsg, "brine")
    if (.not. c_associated(raw)) return
    if (.not. (present(z_enter) .or. present(z_exit) .or. present(thickness) .or. present(found))) then
      call check_status(int(SEAMS_ERR_ARG, c_int), raw, stat, errmsg, "fronts needs one output")
      return
    end if
    enter_ptr = c_null_ptr
    exit_ptr = c_null_ptr
    thick_ptr = c_null_ptr
    found_ptr = c_null_ptr
    if (present(z_enter)) enter_ptr = c_loc(enter_b)
    if (present(z_exit)) exit_ptr = c_loc(exit_b)
    if (present(thickness)) thick_ptr = c_loc(thick_b)
    if (present(found)) found_ptr = c_loc(found_b)
    code = seams_brine_fronts(raw, enter_ptr, exit_ptr, thick_ptr, found_ptr)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    if (present(z_enter)) z_enter = enter_b
    if (present(z_exit)) z_exit = exit_b
    if (present(thickness)) thickness = thick_b
    if (present(found)) found = int(found_b, seams_int)
    call check_status(code, raw, stat, errmsg)
  end procedure brine_fronts
  module procedure brine_rejection
    type(c_ptr) :: raw, ice_ptr, liquid_ptr, k_ptr, defined_ptr
    integer(c_int) :: code
    real(c_double), target :: ice_b, liquid_b, k_b
    integer(c_int), target :: defined_b
    raw = open_raw(self%box, stat, errmsg, "brine")
    if (.not. c_associated(raw)) return
    if (.not. (present(c_ice) .or. present(c_liquid) .or. present(k) .or. present(defined))) then
      call check_status(int(SEAMS_ERR_ARG, c_int), raw, stat, errmsg, "rejection needs one output")
      return
    end if
    ice_ptr = c_null_ptr
    liquid_ptr = c_null_ptr
    k_ptr = c_null_ptr
    defined_ptr = c_null_ptr
    if (present(c_ice)) ice_ptr = c_loc(ice_b)
    if (present(c_liquid)) liquid_ptr = c_loc(liquid_b)
    if (present(k)) k_ptr = c_loc(k_b)
    if (present(defined)) defined_ptr = c_loc(defined_b)
    code = seams_brine_rejection(raw, ice_ptr, liquid_ptr, k_ptr, defined_ptr)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    if (present(c_ice)) c_ice = ice_b
    if (present(c_liquid)) c_liquid = liquid_b
    if (present(k)) k = k_b
    if (present(defined)) defined = int(defined_b, seams_int)
    call check_status(code, raw, stat, errmsg)
  end procedure brine_rejection
  module procedure brine_profile
    type(c_ptr) :: raw
    integer(c_size_t) :: n, cap
    integer(c_int) :: code
    real(seams_real), allocatable, target :: zbuf(:), mbuf(:), fbuf(:)
    raw = open_raw(self%box, stat, errmsg, "brine")
    if (.not. c_associated(raw)) return
    code = seams_brine_profile(raw, c_null_ptr, c_null_ptr, c_null_ptr, 0_c_size_t, n)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    allocate(zbuf(n), mbuf(n), fbuf(n))
    cap = n
    code = seams_brine_profile(raw, c_loc(zbuf), c_loc(mbuf), c_loc(fbuf), cap, n)
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    z_mid = zbuf
    molality = mbuf
    ice_fraction = fbuf
    call check_status(code, raw, stat, errmsg)
  end procedure brine_profile
  module procedure brine_pockets
    type(c_ptr) :: raw
    integer(c_int), target :: code, pockets, waters
    n_pockets = 0_seams_int
    reservoir_waters = 0_seams_int
    raw = open_raw(self%box, stat, errmsg, "brine")
    if (.not. c_associated(raw)) return
    code = seams_brine_pockets(raw, c_loc(pockets), c_loc(waters))
    if (code /= SEAMS_OK) then
      call check_status(code, raw, stat, errmsg)
      return
    end if
    n_pockets = int(pockets, seams_int)
    reservoir_waters = int(waters, seams_int)
    call check_status(code, raw, stat, errmsg)
  end procedure brine_pockets
end submodule seams_methods
