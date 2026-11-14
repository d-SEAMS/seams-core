# Neighbours, rings, and CHILL+ from R.
# The package returns positions, the ring table, and the label column.
# Directed cutoff pairs are counted from those positions with the
# minimum image of the dump's orthorhombic box.

args <- commandArgs(trailingOnly = TRUE)
root <- Sys.getenv("DSEAMS_ROOT", "")
if (!nzchar(root)) stop("DSEAMS_ROOT is unset")
lib <- Sys.getenv("DSEAMS_R_LIB", "")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(dseams))

path <- file.path(root, "input", "traj", "ice1h.lammpstrj")
frame <- read_frame(path, type = -1L)

neighbour_pairs <- function(frame, path, cutoff = 3.5) {
  xyz <- positions(frame)
  header <- readLines(path, n = 8L, warn = FALSE)
  box <- vapply(strsplit(header[6:8], "\\s+"), function(parts) {
    as.numeric(parts[2]) - as.numeric(parts[1])
  }, numeric(1))
  n <- nrow(xyz)
  pairs <- 0L
  cut2 <- cutoff * cutoff
  for (i in seq_len(n)) {
    delta <- sweep(xyz, 2L, xyz[i, ])
    delta <- delta - box[col(delta)] * round(delta / box[col(delta)])
    r2 <- rowSums(delta * delta)
    pairs <- pairs + sum(r2 > 1e-8 & r2 < cut2)
  }
  pairs
}

pairs <- neighbour_pairs(frame, path, 3.5)
ice <- chill_plus(frame, cutoff = 3.5, graph = "cutoff")
counts <- as.integer(attr(ice, "counts"))
rings <- primitive_rings(frame, cutoff = 3.5, graph = "cutoff", max_depth = 6L)
histogram <- integer(7)
if (nrow(rings) > 0L) {
  histogram[rings$size + 1L] <- rings$n
}
cat(
  "pairs", pairs,
  "rings", paste(histogram, collapse = " "),
  "chill", paste(counts, collapse = " "),
  "\n"
)
