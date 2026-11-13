! Copyright (c) d-SEAMS developers
! Eight-cell cubic ice. The counts are the ones quoted in docs/api/fortran.md.
program cubic
  use seams
  implicit none
  integer, parameter :: nrep = 2
  integer, parameter :: n = nrep * nrep * nrep * 8
  real(seams_real) :: xyz(3, n), a, length, fraction, mean, cn
  integer(seams_int), allocatable :: labels(:), counts(:), hist(:)
  real(seams_real), allocatable :: q4(:), q6(:)
  type(seams_frame) :: frame
  type(seams_rings) :: rings
  type(seams_neighbours) :: nb
  integer(seams_int), allocatable :: columns(:)
  real(seams_real) :: fcc(3, 4), basis(3, 2)
  integer :: ii, jj, kk, f, s, placed, i, cubic_n
  character(len=:), allocatable :: version

  fcc = reshape([0.0_seams_real, 0.0_seams_real, 0.0_seams_real, &
    0.5_seams_real, 0.5_seams_real, 0.0_seams_real, &
    0.5_seams_real, 0.0_seams_real, 0.5_seams_real, &
    0.0_seams_real, 0.5_seams_real, 0.5_seams_real], [3, 4])
  basis = reshape([0.0_seams_real, 0.0_seams_real, 0.0_seams_real, &
    0.25_seams_real, 0.25_seams_real, 0.25_seams_real], [3, 2])
  a = 2.75_seams_real * 4.0_seams_real / sqrt(3.0_seams_real)
  placed = 0
  do ii = 0, nrep - 1
    do jj = 0, nrep - 1
      do kk = 0, nrep - 1
        do f = 1, 4
          do s = 1, 2
            placed = placed + 1
            xyz(1, placed) = (ii + fcc(1, f) + basis(1, s)) * a
            xyz(2, placed) = (jj + fcc(2, f) + basis(2, s)) * a
            xyz(3, placed) = (kk + fcc(3, f) + basis(3, s)) * a
          end do
        end do
      end do
    end do
  end do
  length = nrep * a
  call set_num_threads(1_seams_int)
  frame = new_frame()
  call frame%set_ortho(xyz, length, length, length)
  call frame%chill_plus(labels, counts)
  cubic_n = 0
  do i = 1, size(labels)
    if (labels(i) == SEAMS_LABEL_CUBIC) cubic_n = cubic_n + 1
  end do
  call frame%q4q6(q4, q6)
  mean = sum(q6) / size(q6)
  nb = frame%neighbours()
  call nb%columns(columns)
  call nb%close()
  rings = frame%rings(max_depth=6_seams_int)
  call rings%histogram(hist)
  call rings%fraction(6_seams_int, fraction)
  cn = frame%coordination()
  version = version_text()
  if (cubic_n /= 64) stop 1
  if (counts(1) /= 64) stop 2
  if (size(columns) /= 256) stop 3
  if (hist(6) /= 128) stop 4
  if (abs(fraction - 1.0_seams_real) > 1.0e-8_seams_real) stop 5
  if (mean <= 0.3_seams_real) stop 6
  block
    integer(seams_int) :: n_hex, n_cubic, n_water, n_both, n_cages, n_sig
    integer(seams_int) :: n_pairs, n_cation, n_anion, n_subset, largest, n_bonds
    integer(seams_int), allocatable :: cage_labels(:), kinds(:), perfect(:), imperfect(:), sheet(:)
    integer(seams_int), allocatable :: sym_labels(:), fp_counts(:)
    real(seams_real) :: q6_mean, q6_median, percolation, volume, z_enter, z_exit, thickness
    real(seams_real), allocatable :: ql(:), qlbar(:), rr(:), gg(:), zz(:), rho(:)
    integer(seams_int) :: n_ice, found
    type(seams_cages) :: cages
    type(seams_topo) :: fp
    type(seams_brine) :: front
    integer(seams_int), allocatable :: ice(:)
    character(len=:), allocatable :: key
    call frame%ql(6_seams_int, ql, qlbar)
    cages = frame%cages()
    call cages%counts(n_hex, n_cubic, n_water, n_both)
    call cages%close()
    call frame%signature("512", n_cages, n_sig)
    call frame%rdf(rr, gg, volume=volume, rmax=6.0_seams_real, bins=8_seams_int)
    n_bonds = frame%hbonds()
    call frame%density(zz, rho, bins=4_seams_int)
    fp = frame%fingerprint()
    call fp%key(key)
    call fp%rings(fp_counts)
    call fp%close()
    call frame%contact_pairs("1=waterO", n_pairs, n_cation, n_anion)
    call frame%domains("1=waterO", n_subset, largest, percolation)
    call frame%prisms(perfect, imperfect)
    call frame%monolayer(sheet, length * length)
    call frame%site_kinds("1=waterO", kinds)
    call frame%symmetry(sym_labels)
    allocate(ice(frame%count()))
    ice = 1_seams_int
    front = frame%brine(ice)
    call front%fronts(z_enter, z_exit, thickness, found)
    call front%close()
    call frame%q6bar(qlbar, q6_mean, q6_median, n_ice)
    print "(a,a)", "version ", version
    print "(a,i0)", "atoms ", frame%count()
    print "(a,i0)", "cubic ", cubic_n
    print "(a,i0)", "counts_cubic ", counts(1)
    print "(a,f12.8)", "q6_mean ", mean
    print "(a,i0)", "pairs ", size(columns)
    print "(a,i0)", "rings6 ", hist(6)
    print "(a,f12.8)", "fraction6 ", fraction
    print "(a,f12.8)", "coordination ", cn
    print "(a,f12.8)", "ql6_first ", ql(1)
    print "(a,i0,1x,i0,1x,i0,1x,i0)", "cages ", n_hex, n_cubic, n_water, n_both
    print "(a,i0,1x,i0)", "signature512 ", n_cages, n_sig
    print "(a,i0,f12.6)", "rdf_bins ", size(rr), volume
    print "(a,i0)", "hbonds ", n_bonds
    print "(a,i0,f12.6)", "density_bins ", size(zz), rho(1)
    print "(a,i0)", "fingerprint ", len(key)
    print "(a,i0,1x,i0,1x,i0)", "contacts ", n_pairs, n_cation, n_anion
    print "(a,i0,1x,i0,f12.6)", "domain ", n_subset, largest, percolation
    print "(a,i0)", "prisms6 ", perfect(6)
    print "(a,i0)", "monolayer6 ", sheet(6)
    print "(a,i0)", "site0 ", kinds(1)
    print "(a,i0)", "symmetry0 ", sym_labels(1)
    print "(a,i0,f12.6,f12.6)", "brine ", found, z_enter, thickness
    print "(a,i0,f12.8)", "q6bar_ice ", n_ice, q6_mean
  end block
  call rings%close()
  call frame%close()
end program cubic
