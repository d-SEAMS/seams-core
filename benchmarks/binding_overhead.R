# Copyright (c) d-SEAMS developers
# Times the R .Call wrappers in seams_call.c. SEAMS_R_SO is that shared object.

so <- Sys.getenv("SEAMS_R_SO", unset = "")
if (!nzchar(so)) {
  stop("SEAMS_R_SO")
}
dll <- dyn.load(so)
sym <- function(name) getNativeSymbolInfo(name, dll)$address

seams_system_new <- function() .Call(sym("seams_r_system"))
seams_system_read <- function(sys, path) {
  .Call(sym("seams_r_read"), sys, normalizePath(path, mustWork = TRUE), 1L, -1L)
}
seams_system_count <- function(sys) .Call(sym("seams_r_count"), sys)
seams_set_num_threads <- function(n) .Call(sym("seams_r_threads"), as.integer(n))
seams_chill <- function(sys) .Call(sym("seams_r_chill"), sys, 3.5, 0L, 0L, 0L)
seams_steinhardt <- function(sys) .Call(sym("seams_r_steinhardt"), sys, 3.5, 6L, 1L)
seams_neighbours_knn <- function(sys) .Call(sym("seams_r_knn"), sys, 4L, 3.5, 1L)
seams_rings_vertices <- function(sys) .Call(sym("seams_r_vertices"), sys, 3.5, 6L, 0L)

frames <- new.env(parent = emptyenv())
sink_value <- 0
con <- file("stdin")
open(con)
repeat {
  line <- readLines(con, n = 1, warn = FALSE)
  if (length(line) == 0 || !nzchar(line) || startsWith(line, "quit")) break
  parts <- strsplit(line, " ", fixed = TRUE)[[1]]
  threads <- as.integer(parts[[1]])
  op <- parts[[2]]
  path <- parts[[3]]
  Sys.setenv(OMP_NUM_THREADS = threads, PLUMED_NUM_THREADS = threads)
  seams_set_num_threads(threads)
  if (!exists(path, envir = frames, inherits = FALSE)) {
    assign(path, seams_system_new(), envir = frames)
  }
  sys <- get(path, envir = frames, inherits = FALSE)
  if (op != "read" && seams_system_count(sys) == 0) {
    seams_system_read(sys, path)
  }
  t0 <- proc.time()[["elapsed"]]
  if (op == "read") {
    seams_system_read(sys, path)
    sink_value <- sink_value + seams_system_count(sys)
  } else if (op == "knn") {
    sink_value <- sink_value + seams_neighbours_knn(sys)
  } else if (op == "steinhardt") {
    got <- seams_steinhardt(sys)
    sink_value <- sink_value + got[[1]][[1]] + got[[2]][[1]]
  } else if (op == "chill") {
    counts <- seams_chill(sys)
    sink_value <- sink_value + counts[[2]]
  } else if (op == "rings") {
    atoms <- seams_rings_vertices(sys)
    sink_value <- sink_value + length(atoms)
  } else {
    stop(op)
  }
  ms <- (proc.time()[["elapsed"]] - t0) * 1e3
  cat(sprintf("%.9f\n", ms))
  flush(stdout())
}
close(con)
message(sprintf("sink %.6f", sink_value))
