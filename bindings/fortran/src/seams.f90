! Copyright (c) d-SEAMS developers
!
!> Fortran API for d-SEAMS.
!!
!! A program `use`s this module. Handles are derived types. The last
!! finalizer, or an explicit `close`, destroys the C object. Assignment
!! shares the handle. Coordinates are `real(seams_real)` with shape
!! `(3, n)`, x/y/z down the first index, which is the same memory order
!! as the C array. A missing `stat` and a missing `errmsg` make a
!! failed call `error stop` with the library message. Present `stat`
!! receives the `seams.h` code and `errmsg` receives the detail.
module seams
  use seams_kinds, only: seams_int, seams_real
  use seams_support, only: handle_box, handle_share, handle_release
  implicit none
  private

  !> Working precision. It is `real64` and it matches C `double`.
  public :: seams_real
  !> Integer kind for labels, counts, and `stat`. It is `int32`.
  public :: seams_int

  !> The call succeeded.
  integer(seams_int), public, parameter :: SEAMS_OK = 0_seams_int
  !> A required pointer is null.
  integer(seams_int), public, parameter :: SEAMS_ERR_NULL = 1_seams_int
  !> A number is outside the accepted range.
  integer(seams_int), public, parameter :: SEAMS_ERR_RANGE = 2_seams_int
  !> The cell is missing or has no volume.
  integer(seams_int), public, parameter :: SEAMS_ERR_BOX = 3_seams_int
  !> A caller argument is the wrong shape.
  integer(seams_int), public, parameter :: SEAMS_ERR_ARG = 4_seams_int
  !> The handle has no frame yet.
  integer(seams_int), public, parameter :: SEAMS_ERR_STATE = 5_seams_int
  !> An allocation failed.
  integer(seams_int), public, parameter :: SEAMS_ERR_NOMEM = 6_seams_int
  !> The library hit an internal failure.
  integer(seams_int), public, parameter :: SEAMS_ERR_INTERNAL = 7_seams_int
  !> Published CHILL+ rule. Inclusive c_ij edges.
  integer(seams_int), public, parameter :: SEAMS_CHILL_PUBLISHED = 0_seams_int
  !> OVITO CHILL+ rule. Open c_ij edges.
  integer(seams_int), public, parameter :: SEAMS_CHILL_OVITO = 1_seams_int
  !> Pairs inside the cutoff.
  integer(seams_int), public, parameter :: SEAMS_GRAPH_CUTOFF = 0_seams_int
  !> Mutual four-nearest neighbours.
  integer(seams_int), public, parameter :: SEAMS_GRAPH_MUTUAL = 1_seams_int
  !> Cage graph of mutual seeds and union completions.
  integer(seams_int), public, parameter :: SEAMS_GRAPH_SEEDED = 2_seams_int
  !> CHILL cubic ice. Population counts use this index plus one.
  integer(seams_int), public, parameter :: SEAMS_LABEL_CUBIC = 0_seams_int
  !> CHILL hexagonal ice.
  integer(seams_int), public, parameter :: SEAMS_LABEL_HEXAGONAL = 1_seams_int
  !> CHILL liquid water.
  integer(seams_int), public, parameter :: SEAMS_LABEL_WATER = 2_seams_int
  !> CHILL interfacial ice.
  integer(seams_int), public, parameter :: SEAMS_LABEL_INTERFACIAL = 3_seams_int
  !> CHILL clathrate. The population count is `counts(5)`.
  integer(seams_int), public, parameter :: SEAMS_LABEL_CLATHRATE = 4_seams_int
  !> CHILL interfacial clathrate.
  integer(seams_int), public, parameter :: SEAMS_LABEL_INTERCLATHRATE = 5_seams_int
  !> CHILL atom that matched no class.
  integer(seams_int), public, parameter :: SEAMS_LABEL_UNCLASSIFIED = 6_seams_int
  !> CHILL reclassified cubic ice.
  integer(seams_int), public, parameter :: SEAMS_LABEL_RECUBIC = 7_seams_int
  !> CHILL reclassified hexagonal ice.
  integer(seams_int), public, parameter :: SEAMS_LABEL_REHEX = 8_seams_int
  !> Atom in neither cage.
  integer(seams_int), public, parameter :: SEAMS_CAGE_WATER = 0_seams_int
  !> Hexagonal cage only.
  integer(seams_int), public, parameter :: SEAMS_CAGE_HEX = 1_seams_int
  !> Double-diamond cage only.
  integer(seams_int), public, parameter :: SEAMS_CAGE_CUBIC = 2_seams_int
  !> Both cages.
  integer(seams_int), public, parameter :: SEAMS_CAGE_BOTH = 3_seams_int

  !> One frame. The finalizer frees the C frame.
  type, public :: seams_frame
    private
    type(handle_box), pointer :: box => null()
  contains
    procedure, public :: set_ortho => frame_set_ortho
    procedure, public :: set_frame => frame_set_frame
    procedure, public :: read => frame_read
    procedure, public :: count => frame_count
    procedure, public :: close => frame_close
    procedure, public :: chill => frame_chill
    procedure, public :: chill_plus => frame_chill_plus
    procedure, public :: q4q6 => frame_q4q6
    procedure, public :: ql => frame_ql
    procedure, public :: steinhardt => frame_steinhardt
    procedure, public :: q6bar => frame_q6bar
    procedure, public :: neighbours => frame_neighbours
    procedure, public :: knn => frame_knn
    procedure, public :: rings => frame_rings
    procedure, public :: cages => frame_cages
    procedure, public :: signature => frame_signature
    procedure, public :: rdf => frame_rdf
    procedure, public :: running_cn => frame_running_cn
    procedure, public :: coordination => frame_coordination
    procedure, public :: hbonds => frame_hbonds
    procedure, public :: density => frame_density
    procedure, public :: fingerprint => frame_fingerprint
    procedure, public :: contact_pairs => frame_contact_pairs
    procedure, public :: domains => frame_domains
    procedure, public :: prisms => frame_prisms
    procedure, public :: monolayer => frame_monolayer
    procedure, public :: site_kinds => frame_site_kinds
    procedure, public :: symmetry => frame_symmetry
    procedure, public :: brine => frame_brine
    final :: frame_final
  end type seams_frame

  !> Cutoff or neighbour graph. `offsets` and `columns` allocate copies.
  type, public :: seams_neighbours
    private
    type(handle_box), pointer :: box => null()
  contains
    procedure, public :: count => neighbour_count
    procedure, public :: offsets => neighbour_offsets
    procedure, public :: columns => neighbour_columns
    procedure, public :: close => neighbour_close
    final :: neighbour_final
  end type seams_neighbours

  !> Primitive rings. `histogram` is 0-based: index k is rings of size k.
  type, public :: seams_rings
    private
    type(handle_box), pointer :: box => null()
    integer(seams_int) :: max_depth = 0_seams_int
  contains
    procedure, public :: count => ring_count
    procedure, public :: histogram => ring_histogram
    procedure, public :: fraction => ring_fraction
    procedure, public :: sizes => ring_sizes
    procedure, public :: symmetry => ring_symmetry
    procedure, public :: vertices => ring_vertices
    procedure, public :: close => ring_close
    final :: ring_final
  end type seams_rings

  !> Hexagonal and double-diamond membership.
  type, public :: seams_cages
    private
    type(handle_box), pointer :: box => null()
    integer(seams_int) :: natoms = 0_seams_int
  contains
    procedure, public :: counts => cage_counts
    procedure, public :: labels => cage_labels
    procedure, public :: close => cage_close
    final :: cage_final
  end type seams_cages

  !> Topological fingerprint.
  type, public :: seams_topo
    private
    type(handle_box), pointer :: box => null()
  contains
    procedure, public :: key => topo_key
    procedure, public :: atom_key => topo_atom_key
    procedure, public :: rings => topo_rings
    procedure, public :: close => topo_close
    final :: topo_final
  end type seams_topo

  !> Planar ice/brine front.
  type, public :: seams_brine
    private
    type(handle_box), pointer :: box => null()
  contains
    procedure, public :: fronts => brine_fronts
    procedure, public :: rejection => brine_rejection
    procedure, public :: profile => brine_profile
    procedure, public :: pockets => brine_pockets
    procedure, public :: close => brine_close
    final :: brine_final
  end type seams_brine

  !> Empty frame. The caller owns the result.
  public :: new_frame
  !> Library version text.
  public :: version_text
  !> Set the OpenMP width. Zero restores the default.
  public :: set_num_threads
  !> Width the next analysis will request.
  public :: num_threads
  !> Name of a CHILL label code.
  public :: label_name
  !> Name of a status code.
  public :: status_text
  !> Number of frames in a trajectory path.
  public :: path_nframes

  public :: assignment(=)
  interface assignment(=)
    module procedure frame_assign
    module procedure neighbour_assign
    module procedure ring_assign
    module procedure cage_assign
    module procedure topo_assign
    module procedure brine_assign
  end interface

  interface
    module function new_frame(stat, errmsg) result(self)
      type(seams_frame) :: self
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end function

    module function version_text() result(text)
      character(len=:), allocatable :: text
    end function

    module subroutine set_num_threads(n, stat, errmsg)
      integer(seams_int), intent(in) :: n
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    module function num_threads() result(n)
      integer(seams_int) :: n
    end function

    module function label_name(code) result(text)
      integer(seams_int), intent(in) :: code
      character(len=:), allocatable :: text
    end function

    module function status_text(code) result(text)
      integer(seams_int), intent(in) :: code
      character(len=:), allocatable :: text
    end function

    module function path_nframes(path, stat, errmsg) result(nframes)
      character(len=*), intent(in) :: path
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
      integer(seams_int) :: nframes
    end function




    !> `xyz` has shape `(3, n)`. Lengths are angstrom. The origin defaults to zero.
    module subroutine frame_set_ortho(self, xyz, lx, ly, lz, ox, oy, oz, types, ids, stat, errmsg)
      class(seams_frame), intent(inout) :: self
      real(seams_real), intent(in), contiguous :: xyz(:, :)
      real(seams_real), intent(in) :: lx, ly, lz
      real(seams_real), intent(in), optional :: ox, oy, oz
      integer(seams_int), intent(in), optional, contiguous :: types(:), ids(:)
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    !> `box` is ox, oy, oz, lx, ly, lz, xy, xz, yz.
    module subroutine frame_set_frame(self, xyz, box, types, ids, stat, errmsg)
      class(seams_frame), intent(inout) :: self
      real(seams_real), intent(in), contiguous :: xyz(:, :)
      real(seams_real), intent(in) :: box(9)
      integer(seams_int), intent(in), optional, contiguous :: types(:), ids(:)
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    !> `frame` is 1-based. A negative `type_filter` keeps every atom. Zero guesses oxygen.
    module subroutine frame_read(self, path, frame, type_filter, stat, errmsg)
      class(seams_frame), intent(inout) :: self
      character(len=*), intent(in) :: path
      integer(seams_int), intent(in), optional :: frame, type_filter
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    module function frame_count(self) result(n)
      class(seams_frame), intent(in) :: self
      integer(seams_int) :: n
    end function

    module subroutine frame_chill(self, labels, counts, cutoff, typ, graph, rule, stat, errmsg)
      class(seams_frame), intent(inout) :: self
      integer(seams_int), allocatable, intent(out) :: labels(:)
      integer(seams_int), allocatable, intent(out), optional :: counts(:)
      real(seams_real), intent(in), optional :: cutoff
      integer(seams_int), intent(in), optional :: typ, graph, rule
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    module subroutine frame_chill_plus(self, labels, counts, cutoff, typ, graph, rule, stat, errmsg)
      class(seams_frame), intent(inout) :: self
      integer(seams_int), allocatable, intent(out) :: labels(:)
      integer(seams_int), allocatable, intent(out), optional :: counts(:)
      real(seams_real), intent(in), optional :: cutoff
      integer(seams_int), intent(in), optional :: typ, graph, rule
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    module subroutine frame_q4q6(self, q4, q6, cutoff, typ, stat, errmsg)
      class(seams_frame), intent(inout) :: self
      real(seams_real), allocatable, intent(out) :: q4(:), q6(:)
      real(seams_real), intent(in), optional :: cutoff
      integer(seams_int), intent(in), optional :: typ
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    !> `order_l` is 3, 4, 6, 8, or 12. `qlbar` is omitted to skip the average.
    module subroutine frame_ql(self, order_l, ql, qlbar, cutoff, typ, graph, stat, errmsg)
      class(seams_frame), intent(inout) :: self
      integer(seams_int), intent(in) :: order_l
      real(seams_real), allocatable, intent(out) :: ql(:)
      real(seams_real), allocatable, intent(out), optional :: qlbar(:)
      real(seams_real), intent(in), optional :: cutoff
      integer(seams_int), intent(in), optional :: typ, graph
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    module subroutine frame_steinhardt(self, order_l, ql, qlbar, cutoff, typ, stat, errmsg)
      class(seams_frame), intent(inout) :: self
      integer(seams_int), intent(in) :: order_l
      real(seams_real), allocatable, intent(out) :: ql(:)
      real(seams_real), allocatable, intent(out), optional :: qlbar(:)
      real(seams_real), intent(in), optional :: cutoff
      integer(seams_int), intent(in), optional :: typ
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    !> `ice_cut` defaults to 0.36. At least one of the outputs is required.
    module subroutine frame_q6bar(self, q6bar, mean, median, n_ice, cutoff, typ, ice_cut, stat, errmsg)
      class(seams_frame), intent(inout) :: self
      real(seams_real), allocatable, intent(out), optional :: q6bar(:)
      real(seams_real), intent(out), optional :: mean, median
      integer(seams_int), intent(out), optional :: n_ice
      real(seams_real), intent(in), optional :: cutoff, ice_cut
      integer(seams_int), intent(in), optional :: typ
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    module function frame_neighbours(self, cutoff, typ, graph, stat, errmsg) result(list)
      class(seams_frame), intent(inout) :: self
      type(seams_neighbours) :: list
      real(seams_real), intent(in), optional :: cutoff
      integer(seams_int), intent(in), optional :: typ, graph
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end function

    !> `mutual` defaults to true.
    module function frame_knn(self, k, cutoff, typ, mutual, stat, errmsg) result(list)
      class(seams_frame), intent(inout) :: self
      integer(seams_int), intent(in) :: k
      type(seams_neighbours) :: list
      real(seams_real), intent(in), optional :: cutoff
      integer(seams_int), intent(in), optional :: typ
      logical, intent(in), optional :: mutual
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end function

    !> `max_depth` defaults to 6 and must be from 3 through 12.
    module function frame_rings(self, cutoff, typ, graph, max_depth, stat, errmsg) result(rings)
      class(seams_frame), intent(inout) :: self
      type(seams_rings) :: rings
      real(seams_real), intent(in), optional :: cutoff
      integer(seams_int), intent(in), optional :: typ, graph, max_depth
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end function

    !> `graph` defaults to the seeded cage graph. `k` defaults to 4.
    module function frame_cages(self, cutoff, typ, graph, k, complete, stat, errmsg) result(cages)
      class(seams_frame), intent(inout) :: self
      type(seams_cages) :: cages
      real(seams_real), intent(in), optional :: cutoff
      integer(seams_int), intent(in), optional :: typ, graph, k
      logical, intent(in), optional :: complete
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end function

    !> `spec` is a census such as "512" or a name such as "hc".
    module subroutine frame_signature(self, spec, n_cages, n_atoms, cutoff, typ, graph, stat, errmsg)
      class(seams_frame), intent(inout) :: self
      character(len=*), intent(in) :: spec
      integer(seams_int), intent(out) :: n_cages, n_atoms
      real(seams_real), intent(in), optional :: cutoff
      integer(seams_int), intent(in), optional :: typ, graph
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    module subroutine frame_rdf(self, r, g, volume, rmax, bins, type_i, type_j, count, stat, errmsg)
      class(seams_frame), intent(inout) :: self
      real(seams_real), allocatable, intent(out) :: r(:), g(:)
      real(seams_real), intent(out), optional :: volume
      real(seams_real), intent(in), optional :: rmax
      integer(seams_int), intent(in), optional :: bins, type_i, type_j
      integer(seams_int), allocatable, intent(out), optional :: count(:)
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    module subroutine frame_running_cn(self, r, cn, rmax, bins, type_i, type_j, stat, errmsg)
      class(seams_frame), intent(inout) :: self
      real(seams_real), allocatable, intent(out) :: r(:), cn(:)
      real(seams_real), intent(in), optional :: rmax
      integer(seams_int), intent(in), optional :: bins, type_i, type_j
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    module function frame_coordination(self, rmax, bins, type_i, type_j, stat, errmsg) result(cn)
      class(seams_frame), intent(inout) :: self
      real(seams_real), intent(in), optional :: rmax
      integer(seams_int), intent(in), optional :: bins, type_i, type_j
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
      real(seams_real) :: cn
    end function

    module function frame_hbonds(self, cutoff, oxygen_type, hydrogen_type, dist, angle, stat, errmsg) &
        result(n_bonds)
      class(seams_frame), intent(inout) :: self
      real(seams_real), intent(in), optional :: cutoff, dist, angle
      integer(seams_int), intent(in), optional :: oxygen_type, hydrogen_type
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
      integer(seams_int) :: n_bonds
    end function

    !> `axis` is 0, 1, or 2. It defaults to 2.
    module subroutine frame_density(self, z, rho, typ, bins, axis, stat, errmsg)
      class(seams_frame), intent(inout) :: self
      real(seams_real), allocatable, intent(out) :: z(:), rho(:)
      integer(seams_int), intent(in), optional :: typ, bins, axis
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    module function frame_fingerprint(self, cutoff, typ, k, graph, hops, colours, stat, errmsg) &
        result(fp)
      class(seams_frame), intent(inout) :: self
      type(seams_topo) :: fp
      real(seams_real), intent(in), optional :: cutoff
      integer(seams_int), intent(in), optional :: typ, k, graph, hops
      logical, intent(in), optional :: colours
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end function

    module subroutine frame_contact_pairs(self, spec, n_pairs, n_cation, n_anion, stat, errmsg)
      class(seams_frame), intent(inout) :: self
      character(len=*), intent(in) :: spec
      integer(seams_int), intent(out) :: n_pairs, n_cation, n_anion
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    !> `polar` false is the apolar subset.
    module subroutine frame_domains(self, spec, n_subset, largest, percolation, cutoff, polar, stat, errmsg)
      class(seams_frame), intent(inout) :: self
      character(len=*), intent(in) :: spec
      integer(seams_int), intent(out) :: n_subset, largest
      real(seams_real), intent(out) :: percolation
      real(seams_real), intent(in), optional :: cutoff
      logical, intent(in), optional :: polar
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    module subroutine frame_prisms(self, perfect, imperfect, cutoff, typ, stat, errmsg)
      class(seams_frame), intent(inout) :: self
      integer(seams_int), allocatable, intent(out) :: perfect(:), imperfect(:)
      real(seams_real), intent(in), optional :: cutoff
      integer(seams_int), intent(in), optional :: typ
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    module subroutine frame_monolayer(self, counts, sheet_area, cutoff, typ, max_depth, stat, errmsg)
      class(seams_frame), intent(inout) :: self
      integer(seams_int), allocatable, intent(out) :: counts(:)
      real(seams_real), intent(in) :: sheet_area
      real(seams_real), intent(in), optional :: cutoff
      integer(seams_int), intent(in), optional :: typ, max_depth
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    module subroutine frame_site_kinds(self, spec, kinds, stat, errmsg)
      class(seams_frame), intent(inout) :: self
      character(len=*), intent(in) :: spec
      integer(seams_int), allocatable, intent(out) :: kinds(:)
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    !> A non-positive `shell_rmsd` selects the library ceiling.
    module subroutine frame_symmetry(self, labels, shell_rmsd, stat, errmsg)
      class(seams_frame), intent(inout) :: self
      integer(seams_int), allocatable, intent(out) :: labels(:)
      real(seams_real), intent(in), optional :: shell_rmsd
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    !> `ice(i)` nonzero marks atom i as ice. `axis` defaults to 2. `n_bins` 0 selects 20.
    module function frame_brine(self, ice, water_type, ion_type, n_bins, axis, stat, errmsg) result(front)
      class(seams_frame), intent(inout) :: self
      integer(seams_int), intent(in), contiguous :: ice(:)
      type(seams_brine) :: front
      integer(seams_int), intent(in), optional :: water_type, ion_type, n_bins, axis
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end function

    module function neighbour_count(self) result(n)
      class(seams_neighbours), intent(in) :: self
      integer(seams_int) :: n
    end function
    module subroutine neighbour_offsets(self, offsets, stat, errmsg)
      class(seams_neighbours), intent(in) :: self
      integer(seams_int), allocatable, intent(out) :: offsets(:)
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine
    module subroutine neighbour_columns(self, columns, stat, errmsg)
      class(seams_neighbours), intent(in) :: self
      integer(seams_int), allocatable, intent(out) :: columns(:)
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    module function ring_count(self) result(n)
      class(seams_rings), intent(in) :: self
      integer(seams_int) :: n
    end function
    module subroutine ring_histogram(self, counts, stat, errmsg)
      class(seams_rings), intent(in) :: self
      integer(seams_int), allocatable, intent(out) :: counts(:)
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine
    !> `size` 0 counts atoms that sit on any ring.
    module subroutine ring_fraction(self, size, fraction, stat, errmsg)
      class(seams_rings), intent(in) :: self
      integer(seams_int), intent(in) :: size
      real(seams_real), intent(out) :: fraction
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine
    module subroutine ring_sizes(self, sizes, stat, errmsg)
      class(seams_rings), intent(in) :: self
      integer(seams_int), allocatable, intent(out) :: sizes(:)
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine
    module subroutine ring_symmetry(self, codes, stat, errmsg)
      class(seams_rings), intent(in) :: self
      integer(seams_int), allocatable, intent(out) :: codes(:)
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine
    !> `index` is 1-based.
    module subroutine ring_vertices(self, index, atoms, stat, errmsg)
      class(seams_rings), intent(in) :: self
      integer(seams_int), intent(in) :: index
      integer(seams_int), allocatable, intent(out) :: atoms(:)
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    module subroutine cage_counts(self, n_hex, n_cubic, n_water, n_both, stat, errmsg)
      class(seams_cages), intent(in) :: self
      integer(seams_int), intent(out) :: n_hex, n_cubic, n_water, n_both
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine
    module subroutine cage_labels(self, labels, stat, errmsg)
      class(seams_cages), intent(in) :: self
      integer(seams_int), allocatable, intent(out) :: labels(:)
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    module subroutine topo_key(self, key, stat, errmsg)
      class(seams_topo), intent(in) :: self
      character(len=:), allocatable, intent(out) :: key
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine
    !> `index` is 1-based.
    module subroutine topo_atom_key(self, index, key, stat, errmsg)
      class(seams_topo), intent(in) :: self
      integer(seams_int), intent(in) :: index
      character(len=:), allocatable, intent(out) :: key
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine
    module subroutine topo_rings(self, counts, stat, errmsg)
      class(seams_topo), intent(in) :: self
      integer(seams_int), allocatable, intent(out) :: counts(:)
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine

    module subroutine brine_fronts(self, z_enter, z_exit, thickness, found, stat, errmsg)
      class(seams_brine), intent(in) :: self
      real(seams_real), intent(out), optional :: z_enter, z_exit, thickness
      integer(seams_int), intent(out), optional :: found
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine
    module subroutine brine_rejection(self, c_ice, c_liquid, k, defined, stat, errmsg)
      class(seams_brine), intent(in) :: self
      real(seams_real), intent(out), optional :: c_ice, c_liquid, k
      integer(seams_int), intent(out), optional :: defined
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine
    module subroutine brine_profile(self, z_mid, molality, ice_fraction, stat, errmsg)
      class(seams_brine), intent(in) :: self
      real(seams_real), allocatable, intent(out) :: z_mid(:), molality(:), ice_fraction(:)
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine
    module subroutine brine_pockets(self, n_pockets, reservoir_waters, stat, errmsg)
      class(seams_brine), intent(in) :: self
      integer(seams_int), intent(out) :: n_pockets, reservoir_waters
      integer(seams_int), intent(out), optional :: stat
      character(len=:), allocatable, intent(out), optional :: errmsg
    end subroutine
  end interface

contains

  !> Share the C frame. The last finalizer destroys it.
  subroutine frame_assign(lhs, rhs)
    type(seams_frame), intent(inout) :: lhs
    type(seams_frame), intent(in) :: rhs
    call handle_share(lhs%box, rhs%box)
  end subroutine frame_assign

  !> Release this owner's reference to the frame.
  subroutine frame_final(self)
    type(seams_frame), intent(inout) :: self
    call handle_release(self%box)
  end subroutine frame_final

  !> Release this owner's reference to the frame.
  subroutine frame_close(self)
    class(seams_frame), intent(inout) :: self
    call handle_release(self%box)
  end subroutine frame_close

  !> Share the neighbour list.
  subroutine neighbour_assign(lhs, rhs)
    type(seams_neighbours), intent(inout) :: lhs
    type(seams_neighbours), intent(in) :: rhs
    call handle_share(lhs%box, rhs%box)
  end subroutine neighbour_assign

  !> Release this owner's reference to the neighbour list.
  subroutine neighbour_final(self)
    type(seams_neighbours), intent(inout) :: self
    call handle_release(self%box)
  end subroutine neighbour_final

  !> Release this owner's reference to the neighbour list.
  subroutine neighbour_close(self)
    class(seams_neighbours), intent(inout) :: self
    call handle_release(self%box)
  end subroutine neighbour_close

  !> Share the ring set and its maximum depth.
  subroutine ring_assign(lhs, rhs)
    type(seams_rings), intent(inout) :: lhs
    type(seams_rings), intent(in) :: rhs
    call handle_share(lhs%box, rhs%box)
    if (associated(rhs%box)) then
      lhs%max_depth = rhs%max_depth
    else
      lhs%max_depth = 0_seams_int
    end if
  end subroutine ring_assign

  !> Release this owner's reference to the ring set.
  subroutine ring_final(self)
    type(seams_rings), intent(inout) :: self
    call handle_release(self%box)
  end subroutine ring_final

  !> Release this owner's reference to the ring set.
  subroutine ring_close(self)
    class(seams_rings), intent(inout) :: self
    call handle_release(self%box)
  end subroutine ring_close

  !> Share the cage set and the atom count it was built for.
  subroutine cage_assign(lhs, rhs)
    type(seams_cages), intent(inout) :: lhs
    type(seams_cages), intent(in) :: rhs
    call handle_share(lhs%box, rhs%box)
    if (associated(rhs%box)) then
      lhs%natoms = rhs%natoms
    else
      lhs%natoms = 0_seams_int
    end if
  end subroutine cage_assign

  !> Release this owner's reference to the cage set.
  subroutine cage_final(self)
    type(seams_cages), intent(inout) :: self
    call handle_release(self%box)
  end subroutine cage_final

  !> Release this owner's reference to the cage set.
  subroutine cage_close(self)
    class(seams_cages), intent(inout) :: self
    call handle_release(self%box)
  end subroutine cage_close

  !> Share the fingerprint.
  subroutine topo_assign(lhs, rhs)
    type(seams_topo), intent(inout) :: lhs
    type(seams_topo), intent(in) :: rhs
    call handle_share(lhs%box, rhs%box)
  end subroutine topo_assign

  !> Release this owner's reference to the fingerprint.
  subroutine topo_final(self)
    type(seams_topo), intent(inout) :: self
    call handle_release(self%box)
  end subroutine topo_final

  !> Release this owner's reference to the fingerprint.
  subroutine topo_close(self)
    class(seams_topo), intent(inout) :: self
    call handle_release(self%box)
  end subroutine topo_close

  !> Share the brine front.
  subroutine brine_assign(lhs, rhs)
    type(seams_brine), intent(inout) :: lhs
    type(seams_brine), intent(in) :: rhs
    call handle_share(lhs%box, rhs%box)
  end subroutine brine_assign

  !> Release this owner's reference to the brine front.
  subroutine brine_final(self)
    type(seams_brine), intent(inout) :: self
    call handle_release(self%box)
  end subroutine brine_final

  !> Release this owner's reference to the brine front.
  subroutine brine_close(self)
    class(seams_brine), intent(inout) :: self
    call handle_release(self%box)
  end subroutine brine_close
end module seams
