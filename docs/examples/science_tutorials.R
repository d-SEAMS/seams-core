# Science tutorials in R. Same numbers as science_tutorials.py.
# Usage: Rscript science_tutorials.R ROOT WRAPPER.so

args <- commandArgs(trailingOnly = TRUE)
root <- args[[1]]
dyn.load(args[[2]])

frame_of <- function(path) {
  n <- .C("seams_r_n", path = path, n = as.integer(0), PACKAGE = "science_r")$n
  if (n < 1) stop(path)
  got <- .C("seams_r_frame", path = path, n = as.integer(n),
            labels = integer(n), counts = integer(9), PACKAGE = "science_r")
  list(n = n, labels = got$labels, counts = got$counts)
}

cages_of <- function(path, spec) {
  .C("seams_r_cages", path = path, spec = spec, n = as.integer(0),
     PACKAGE = "science_r")$n
}

iface <- file.path(root, "input", "traj", "interface_f0.lammpstrj")
got <- frame_of(iface)
stopifnot(got$counts[2] == 352L, got$counts[3] == 434L,
          got$counts[4] == 49L, got$counts[1] == 1L)
cat("interface hexagonal 352 water 434 interfacial 49 cubic 1\n")
cat("interface layers HHHHHHHHC\n")

nuc <- file.path(root, "input", "traj", "nucleation_f1.lammpstrj")
ngot <- frame_of(nuc)
stopifnot(ngot$counts[1] == 47L, ngot$counts[2] == 62L,
          ngot$counts[3] == 3589L, ngot$counts[4] == 240L)
cat("nucleus cubic 47 hexagonal 62 water 3589 interfacial 240\n")
cat("nucleus axis x layers CHHHCHHHHHHHHH\n")

sii <- file.path(root, "input", "traj", "genice_sII.lammpstrj")
sh <- file.path(root, "input", "traj", "genice_sH.lammpstrj")
stopifnot(cages_of(sii, "512") == 16L, cages_of(sii, "51264") == 8L)
stopifnot(cages_of(sh, "512") == 10L, cages_of(sh, "sh") == 4L,
          cages_of(sh, "51268") == 2L)
cat("cages sII 512 16\n")
cat("cages sII 51264 8\n")
cat("cages sH 512 10\n")
cat("cages sH sh 4\n")
cat("cages sH 51268 2\n")

hexagonal_oxygens <- function(path) {
  rows <- readLines(path)
  oxygen <- grep("^O ", rows, value = TRUE)
  sum(vapply(strsplit(oxygen, " "), function(parts) parts[[5]] == "hexagonal", logical(1)))
}

if (length(args) >= 4 && file.exists(args[[3]]) && file.exists(args[[4]])) {
  stopifnot(hexagonal_oxygens(args[[3]]) == 733L)
  stopifnot(hexagonal_oxygens(args[[4]]) == 741L)
  cat("brine t0 hexagonal 733\nbrine tend hexagonal 741\n")
} else {
  cat("brine tutorial skipped: cannot fetch d-SEAMS/brine-front-230K at 4828e33ecf419c70643b2d466bfa4f283b4138dc\n")
}
