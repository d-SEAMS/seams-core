The cutoff cell window is timed with `steady_clock` when OpenMP is off, so the macOS R package and the wheels build that loop without `omp_get_wtime`.
