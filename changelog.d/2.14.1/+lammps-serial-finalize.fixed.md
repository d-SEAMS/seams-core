The serial LAMMPS hook calls `lammps_mpi_finalize` after `lammps_close`, which frees the ReadData and Run commands the stub would otherwise leave allocated. The MPI build still calls `MPI_Finalize`.
