! Copyright (c) d-SEAMS developers
! The 64-atom cubic cell and the trajectory counts match test_seams.f90.
program test_high
  use seams
  implicit none
  integer, parameter :: nrep = 2
  integer, parameter :: n = nrep * nrep * nrep * 8
  real(seams_real) :: xyz(3, n), a, length, fraction, cn
  real(seams_real) :: box(9)
  real(seams_real), allocatable :: q4(:), q6(:)
  integer(seams_int), allocatable :: labels(:), counts(:), hist(:), atoms(:)
  integer(seams_int) :: i, cubic, stat
  integer(seams_int) :: n_hex, n_cubic, n_water, n_both, n_bonds
  type(seams_frame) :: frame
  type(seams_neighbours) :: nb, knn
  type(seams_rings) :: rings
  type(seams_cages) :: cages
  type(seams_topo) :: fp
  real(seams_real) :: fcc(3, 4), basis(3, 2)
  integer :: ii, jj, kk, f, s, placed
  integer(seams_int), allocatable :: columns(:), kinds(:)
  real(seams_real), allocatable :: rr(:), gg(:), zz(:), rho(:)
  integer(seams_int) :: n_cages, n_atoms_sig, n_pairs, n_cation, n_anion
  integer(seams_int) :: n_subset, largest
  real(seams_real) :: percolation, volume
  integer(seams_int), allocatable :: perfect(:), imperfect(:), sheet(:), fp_counts(:)
  character(len=512) :: root
  integer :: nroot
  character(len=:), allocatable :: version, key, err

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
  box = [0.0_seams_real, 0.0_seams_real, 0.0_seams_real, length, length, length, &
    0.0_seams_real, 0.0_seams_real, 0.0_seams_real]

  version = version_text()
  if (len(version) == 0) stop 1
  call set_num_threads(1_seams_int, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 2
  if (num_threads() /= 1) stop 3
  frame = new_frame(stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 4
  call frame%set_ortho(xyz, length, length, length, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 5
  call frame%set_frame(xyz, box, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 6
  if (frame%count() /= n) stop 7
  call frame%chill_plus(labels, counts, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 8
  cubic = 0
  do i = 1, n
    if (labels(i) == SEAMS_LABEL_CUBIC) cubic = cubic + 1
  end do
  if (cubic /= n) stop 9
  call frame%chill_plus(labels, counts, graph=SEAMS_GRAPH_MUTUAL, rule=SEAMS_CHILL_OVITO, &
    stat=stat, errmsg=err)
  if (stat /= SEAMS_OK .or. counts(1) /= n) stop 10
  call frame%q4q6(q4, q6, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 11
  if (sum(q6) / n <= 0.3_seams_real) stop 12
  nb = frame%neighbours(stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 13
  call nb%columns(columns, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK .or. size(columns) <= 0) stop 14
  call nb%close()
  rings = frame%rings(max_depth=6_seams_int, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 15
  call rings%histogram(hist, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK .or. hist(6) <= 0) stop 16
  call rings%fraction(6_seams_int, fraction, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK .or. fraction <= 0.5_seams_real) stop 17
  if (rings%count() <= 0) stop 18
  call rings%vertices(1_seams_int, atoms, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK .or. size(atoms) <= 0) stop 19
  call rings%close()
  if (command_argument_count() < 1) stop 30
  call get_command_argument(1, root, nroot)
  call frame%chill(labels, counts, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 40
  call frame%steinhardt(6_seams_int, q4, q6, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 41
  cages = frame%cages(stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 42
  call cages%counts(n_hex, n_cubic, n_water, n_both, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 43
  call cages%close()
  call frame%signature("512", n_cages, n_atoms_sig, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 44
  knn = frame%knn(4_seams_int, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 45
  call knn%close()
  call frame%rdf(rr, gg, volume=volume, rmax=6.0_seams_real, bins=8_seams_int, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK .or. size(rr) <= 0) stop 46
  call frame%running_cn(rr, gg, rmax=6.0_seams_real, bins=8_seams_int, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK .or. size(rr) <= 0) stop 47
  cn = frame%coordination(rmax=3.5_seams_real, bins=8_seams_int, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 48
  n_bonds = frame%hbonds(stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 49
  call frame%density(zz, rho, bins=4_seams_int, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK .or. size(zz) <= 0) stop 50
  fp = frame%fingerprint(stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 51
  call fp%key(key, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK .or. len(key) <= 0) stop 52
  call fp%rings(fp_counts, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 53
  call fp%close()
  call frame%contact_pairs("1=cationHead", n_pairs, n_cation, n_anion, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 54
  call frame%domains("1=waterO", n_subset, largest, percolation, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 55
  call frame%prisms(perfect, imperfect, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 56
  call frame%monolayer(sheet, 100.0_seams_real, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 57
  call frame%site_kinds("1=waterO", kinds, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) stop 60
  call frame%read(trim(root) // "/input/traj/ice1h.lammpstrj", stat=stat, errmsg=err)
  if (stat /= SEAMS_OK .or. frame%count() /= 128) stop 61
  call frame%read(trim(root) // "/input/traj/genice_sI.lammpstrj", type_filter=1_seams_int, &
    stat=stat, errmsg=err)
  if (stat /= SEAMS_OK .or. frame%count() /= 46) stop 62
  call frame%read(trim(root) // "/input/xyz/configuration.xyz", stat=stat, errmsg=err)
  if (stat /= SEAMS_OK .or. frame%count() <= 0) stop 63
  call frame%read(trim(root) // "/input/con/tiny_multi_cuh2.con", stat=stat, errmsg=err)
  if (stat /= SEAMS_OK .or. frame%count() <= 0) stop 64
  call frame%read(trim(root) // "/input/traj/genice_sH.lammpstrj", type_filter=1_seams_int, &
    stat=stat, errmsg=err)
  if (stat /= SEAMS_OK .or. frame%count() <= 0) stop 65
  call frame%close()
  print *, "fortran high ok"
end program test_high
