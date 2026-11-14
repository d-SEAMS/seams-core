! Copyright (c) d-SEAMS developers
program test_seams
  use, intrinsic :: iso_c_binding
  use seams_c
  implicit none
  integer, parameter :: nrep = 2
  integer, parameter :: n = nrep * nrep * nrep * 8
  real(c_double) :: xyz(3 * n), q4(n), q6(n), length, a, fraction
  integer(c_int) :: labels(n), counts(9), code, i, cubic
  integer(c_int) :: n_offsets, nnz
  real(c_double) :: box9(9)
  type(c_ptr) :: sys, nb, rings, offsets, cols, ver
  integer(c_int) :: hist(0:6)
  real(c_double) :: fcc(3, 4), basis(3, 2)
  integer :: ii, jj, kk, f, s, cursor
  character(len=512) :: root
  integer :: nroot
  integer(c_int) :: n_bonds, n_hex, n_cubic, n_water, n_both
  integer(c_int) :: perfect(7), imperfect(7), sheet(7), kinds(64), fp_counts(13), atoms(12)
  real(c_double), target :: qr(64), qbar(64), rr(64), gg(64), zz(64), rho(64)
  real(c_double), target :: volume, cn, percolation
  integer(c_size_t) :: n_out, n_cages, n_atoms_sig
  type(c_ptr) :: cages, knn, fp

  fcc = reshape([0d0, 0d0, 0d0, 0.5d0, 0.5d0, 0d0, 0.5d0, 0d0, 0.5d0, 0d0, 0.5d0, 0.5d0], [3, 4])
  basis = reshape([0d0, 0d0, 0d0, 0.25d0, 0.25d0, 0.25d0], [3, 2])
  a = 2.75d0 * 4d0 / sqrt(3d0)
  cursor = 0
  do ii = 0, nrep - 1
    do jj = 0, nrep - 1
      do kk = 0, nrep - 1
        do f = 1, 4
          do s = 1, 2
            cursor = cursor + 1
            xyz(3 * cursor - 2) = (ii + fcc(1, f) + basis(1, s)) * a
            xyz(3 * cursor - 1) = (jj + fcc(2, f) + basis(2, s)) * a
            xyz(3 * cursor) = (kk + fcc(3, f) + basis(3, s)) * a
          end do
        end do
      end do
    end do
  end do
  length = nrep * a

  ver = seams_version_string()
  if (.not. c_associated(ver)) stop 14
  ver = seams_last_error()
  if (seams_set_num_threads(1) /= SEAMS_OK) stop 1
  if (seams_get_num_threads() /= 1) stop 2
  sys = seams_system_create()
  if (.not. c_associated(sys)) stop 3
  code = seams_system_set_ortho(sys, xyz, n, 0d0, 0d0, 0d0, length, length, length, &
                                c_null_ptr, c_null_ptr)
  if (code /= SEAMS_OK) stop 4
  box9 = [0d0, 0d0, 0d0, length, length, length, 0d0, 0d0, 0d0]
  code = seams_system_set_frame(sys, xyz, n, box9, c_null_ptr, c_null_ptr)
  if (code /= SEAMS_OK) stop 15
  code = seams_chill_plus(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, counts)
  if (code /= SEAMS_OK) stop 5
  cubic = 0
  do i = 1, n
    if (labels(i) == SEAMS_LABEL_CUBIC) cubic = cubic + 1
  end do
  if (cubic /= n) stop 6
  code = seams_chill_plus(sys, 3.5d0, 0, SEAMS_GRAPH_MUTUAL, SEAMS_CHILL_OVITO, labels, counts)
  if (code /= SEAMS_OK .or. counts(1) /= n) stop 7
  code = seams_q4q6(sys, 3.5d0, 0, q4, q6)
  if (code /= SEAMS_OK) stop 8
  if (sum(q6) / n <= 0.3d0) stop 9
  print '(a,i0)', 'fortran call cubic ', cubic
  print '(a,f14.12)', 'fortran call q6 ', sum(q6) / n
  code = seams_neighbours_build(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, nb)
  if (code /= SEAMS_OK) stop 10
  code = seams_neighbours_csr(nb, offsets, n_offsets, cols, nnz)
  if (code /= SEAMS_OK .or. nnz <= 0) stop 16
  call seams_neighbours_destroy(nb)
  code = seams_rings_find(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, 6, rings)
  if (code /= SEAMS_OK) stop 11
  code = seams_rings_histogram(rings, hist, 7)
  if (code /= SEAMS_OK .or. hist(6) <= 0) stop 12
  code = seams_rings_fraction(rings, 6, fraction)
  if (code /= SEAMS_OK .or. fraction <= 0.5d0) stop 13
  call seams_rings_destroy(rings)
  if (command_argument_count() < 1) stop 30
  call get_command_argument(1, root, nroot)
  code = seams_chill(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, &
                     int(n, c_size_t), counts)
  if (code /= SEAMS_OK) stop 40
  code = seams_steinhardt(sys, 3.5d0, 0, 6, qr, int(n, c_size_t), qbar, int(n, c_size_t))
  if (code /= SEAMS_OK) stop 41
  if (.not. c_associated(seams_handle_error(sys))) stop 41
  code = seams_cages(sys, 3.5d0, 0, SEAMS_GRAPH_SEEDED, 4, 0, cages)
  if (code /= SEAMS_OK) stop 42
  code = seams_cages_counts(cages, n_hex, n_cubic, n_water, n_both)
  if (code /= SEAMS_OK) stop 43
  print '(a,i0,1x,i0,1x,i0,1x,i0)', 'fortran call cages ', n_hex, n_cubic, n_water, n_both
  call seams_cages_destroy(cages)
  code = signature_count(sys, "512", n_cages, n_atoms_sig)
  if (code /= SEAMS_OK) stop 44
  print '(a,i0,1x,i0)', 'fortran call sig ', n_cages, n_atoms_sig
  code = seams_neighbours_knn(sys, 4, 3.5d0, 0, 1, knn)
  if (code /= SEAMS_OK) stop 45
  call seams_neighbours_destroy(knn)
  code = seams_rdf(sys, 6d0, 8, 0, 0, c_loc(rr), c_loc(gg), c_null_ptr, 64_c_size_t, n_out, &
                   c_loc(volume))
  if (code /= SEAMS_OK .or. n_out <= 0) stop 46
  code = seams_running_cn(sys, 6d0, 8, 0, 0, c_loc(rr), c_loc(gg), 64_c_size_t, n_out)
  if (code /= SEAMS_OK .or. n_out <= 0) stop 47
  code = seams_coordination(sys, 3.5d0, 8, 0, 0, cn)
  if (code /= SEAMS_OK) stop 48
  print '(a,f14.12)', 'fortran call cn ', cn
  code = seams_hbonds(sys, 3.5d0, 1, 2, 2.42d0, 30d0, n_bonds)
  if (code /= SEAMS_OK) stop 49
  print '(a,i0)', 'fortran call hbond ', n_bonds
  code = seams_density_axis(sys, 0, 4, 2, c_loc(zz), c_loc(rho), 64_c_size_t, n_out)
  if (code /= SEAMS_OK .or. n_out <= 0) stop 50
  code = seams_fingerprint(sys, 3.5d0, 0, 4, SEAMS_GRAPH_CUTOFF, 2, 0, fp)
  if (code /= SEAMS_OK) stop 51
  code = seams_fingerprint_key(fp, c_null_ptr, 0_c_size_t, n_out)
  if (code /= SEAMS_OK .or. n_out <= 0) stop 52
  code = seams_fingerprint_rings(fp, fp_counts, 13_c_size_t)
  if (code /= SEAMS_OK) stop 53
  call seams_fingerprint_destroy(fp)
  code = pair_count(sys, "1=cationHead", n_cages, n_hex, n_cubic)
  if (code /= SEAMS_OK) stop 54
  code = domain_size(sys, "1=waterO", n_hex, n_cubic, percolation)
  if (code /= SEAMS_OK) stop 55
  code = seams_prisms(sys, 3.5d0, 0, perfect, imperfect, 7_c_size_t)
  if (code /= SEAMS_OK) stop 56
  code = seams_monolayer(sys, 3.5d0, 0, 6, 100d0, sheet, 7_c_size_t)
  if (code /= SEAMS_OK) stop 57
  print '(a,i0)', 'fortran call sheet ', sheet(7)
  code = seams_rings_find(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, 6, rings)
  if (code /= SEAMS_OK) stop 58
  if (seams_rings_count(rings) <= 0) stop 58
  code = seams_rings_vertices(rings, 0, atoms, 12_c_size_t, n_out)
  if (code /= SEAMS_OK .or. n_out <= 0) stop 59
  print '(a,i0,1x,i0)', 'fortran call vert ', n_out, atoms(1)
  call seams_rings_destroy(rings)
  code = site_kinds(sys, "1=waterO", kinds, int(n, c_size_t))
  if (code /= SEAMS_OK) stop 60
  code = read_frame(sys, trim(root) // "/input/traj/ice1h.lammpstrj", 1, -1)
  if (code /= SEAMS_OK) stop 61
  if (seams_system_count(sys) /= 128) stop 61
  code = read_frame(sys, trim(root) // "/input/traj/genice_sI.lammpstrj", 1, 1)
  if (code /= SEAMS_OK) stop 62
  if (seams_system_count(sys) /= 46) stop 62
  code = read_frame(sys, trim(root) // "/input/xyz/configuration.xyz", 1, -1)
  if (code /= SEAMS_OK) stop 63
  if (seams_system_count(sys) <= 0) stop 63
  code = read_frame(sys, trim(root) // "/input/con/tiny_multi_cuh2.con", 1, -1)
  if (code /= SEAMS_OK) stop 64
  if (seams_system_count(sys) <= 0) stop 64
  code = read_frame(sys, trim(root) // "/input/traj/genice_sH.lammpstrj", 1, 1)
  if (code /= SEAMS_OK) stop 65
  if (seams_system_count(sys) <= 0) stop 65
  call seams_system_destroy(sys)
  print '(a)', 'fortran call read 128 46'
  print *, "fortran ok"
end program test_seams
