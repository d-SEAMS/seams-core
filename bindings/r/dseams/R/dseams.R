#' Read one trajectory frame.
#'
#' `type = 0` guesses oxygen (LAMMPS type 2, otherwise type 1). A positive
#' type keeps that species. A negative type keeps every atom.
#'
#' @param path Path to a LAMMPS dump, or an XYZ file.
#' @param frame 1-based frame index.
#' @param type LAMMPS type filter.
#' @return An object of class `dseams_frame`.
#' @export
read_frame <- function(path, frame = 1L, type = 0L) {
  path <- normalizePath(path, mustWork = TRUE)
  ptr <- read_frame_(path, as.integer(frame), as.integer(type))
  structure(
    list(ptr = ptr, path = path, frame = as.integer(frame), type = as.integer(type)),
    class = "dseams_frame"
  )
}

#' Number of LAMMPS timesteps in a dump. An XYZ file counts as one.
#' @param path Trajectory path.
#' @export
nframes <- function(path) {
  nframes_(normalizePath(path, mustWork = TRUE))
}

#' @export
print.dseams_frame <- function(x, ...) {
  cat("dseams frame", x$frame, "nop", nop_(x$ptr), "\n", sep = " ")
  invisible(x)
}

#' Atom positions, one row per atom.
#' @param frame A `dseams_frame`.
#' @export
positions <- function(frame) {
  UseMethod("positions")
}

#' @export
positions.dseams_frame <- function(frame) {
  n <- nop_(frame$ptr)
  matrix(positions_(frame$ptr), nrow = n, ncol = 3, byrow = TRUE,
         dimnames = list(NULL, c("x", "y", "z")))
}

#' CHILL+ ice labels.
#'
#' `graph = "cutoff"` (the default) bonds every pair closer than `cutoff`,
#' as `seams chill-plus`, pydseamslib, and Julia do by default.
#' `graph = "mutual"` keeps the mutual four nearest neighbours inside
#' `cutoff + 1.5`.
#'
#' `rule = "chill+"` is the published rule: the cutoff shell must hold
#' exactly four neighbours, q3 is averaged over the four nearest, and the
#' Table 1 interfacial-neighbour test is on. `rule = "ovito-chill+"`
#' averages q3 over the whole shell and skips that test.
#'
#' @param frame A `dseams_frame`.
#' @param cutoff Neighbour cutoff in angstrom. For the mutual graph this is
#'   only the cell-size hint (`cutoff + 1.5`).
#' @param type LAMMPS type. 0 uses the first atom's type.
#' @param graph `"cutoff"` or `"mutual"`.
#' @param rule `"chill+"` or `"ovito-chill+"`. A YAML or Lua snippet with
#'   a `rule` field is accepted.
#' @return A data.frame with columns `atom`, `type`, `label`, `name`.
#'   The `counts` attribute is the nine-label histogram.
#' @export
chill_plus <- function(frame, cutoff = 3.5, type = 0L, graph = "cutoff",
                       rule = "chill+") {
  UseMethod("chill_plus")
}

#' @export
chill_plus.dseams_frame <- function(frame, cutoff = 3.5, type = 0L,
                                    graph = "cutoff", rule = "chill+") {
  out <- chill_plus_(frame$ptr, cutoff, as.integer(type), graph, rule)
  class(out) <- c("dseams_ice", "data.frame")
  out
}

#' Hexagonal-cage and double-diamond membership.
#'
#' @param frame A `dseams_frame`.
#' @param cutoff Neighbour cutoff in angstrom.
#' @param type LAMMPS type. 0 uses the first atom's type.
#' @param k Neighbours kept by the k-nearest graphs.
#' @param graph `"seeded"`, `"mutual"`, `"cutoff"`, or `"knn-union"`.
#' @return A data.frame with columns `atom`, `cage`, `name`.
#'   Attributes `n_ih`, `n_ic`, `n_both`, `n_water` are disjoint and match
#'   `seams cages` (hexagonal only, cubic only, both, water). Per-atom
#'   `cage` is 0 water, 1 hexagonal, 2 cubic, 3 both.
#' @export
cages <- function(frame, cutoff = 3.5, type = 0L, k = 4L, graph = "seeded") {
  UseMethod("cages")
}

#' @export
cages.dseams_frame <- function(frame, cutoff = 3.5, type = 0L, k = 4L,
                               graph = "seeded") {
  out <- cages_(frame$ptr, cutoff, as.integer(type), as.integer(k), graph)
  class(out) <- c("dseams_cages", "data.frame")
  out
}

#' Primitive ring census.
#'
#' @param frame A `dseams_frame`.
#' @param cutoff Neighbour cutoff in angstrom.
#' @param type LAMMPS type.
#' @param graph `"mutual"` or `"cutoff"`.
#' @param max_depth Largest ring retained.
#' @return A data.frame with columns `size` and `n`.
#' @export
primitive_rings <- function(frame, cutoff = 3.5, type = 0L, graph = "mutual",
                            max_depth = 6L) {
  UseMethod("primitive_rings")
}

#' Closed cages for one ring census.
#'
#' @param frame A `dseams_frame`.
#' @param spec Ring census such as `"512"` or `"5:12,6:2"`.
#' @return An integer vector `c(n_cages, n_atoms)`.
#' @export
signature_cages <- function(frame, spec = "512", cutoff = 3.5, type = 0L,
                            graph = "cutoff") {
  UseMethod("signature_cages")
}

#' @export
signature_cages.dseams_frame <- function(frame, spec = "512", cutoff = 3.5,
                                         type = 0L, graph = "cutoff") {
  signature_cages_(frame$ptr, cutoff, as.integer(type), graph, spec)
}

#' Jump-rotor tau_90 between two frames.
#'
#' Returns `dt` when any water's H-H vector rotates by at least 90 degrees,
#' and `NaN` otherwise.
#' @export
jump_tau90 <- function(first, second, dt = 1, oxygen_type = 1L,
                       hydrogen_type = 2L) {
  jump_tau90_(first$ptr, second$ptr, dt, as.integer(oxygen_type),
              as.integer(hydrogen_type))
}

#' Shell symmetry label for every atom.
#'
#' @return Integer codes: 0 ice Ih, 1 ice Ic, 2 stacking-disordered, -1 otherwise.
#' @export
symmetry <- function(frame, shell_rmsd = 0.12) {
  UseMethod("symmetry")
}

#' @export
symmetry.dseams_frame <- function(frame, shell_rmsd = 0.12) {
  symmetry_(frame$ptr, shell_rmsd)
}

#' @export
primitive_rings.dseams_frame <- function(frame, cutoff = 3.5, type = 0L,
                                         graph = "mutual", max_depth = 6L) {
  out <- rings_(frame$ptr, cutoff, as.integer(type), graph, as.integer(max_depth))
  class(out) <- c("dseams_rings", "data.frame")
  out
}

#' Steinhardt order parameters.
#'
#' `ql` is the local parameter. `qlbar` is the Lechner-Dellago neighbour average.
#'
#' @param frame A `dseams_frame`.
#' @param order_l Spherical-harmonic degree. 3, 4, 6, 8, and 12 are computed.
#'   Degrees 4 and 6 on an orthorhombic or restricted-triclinic cell use the
#'   solid-harmonic kernel.
#' @param cutoff Neighbour cutoff in angstrom. An orthorhombic or
#'   restricted-triclinic cell uses the cutoff shell. Any other cell keeps
#'   the mutual four nearest inside `cutoff + 1.5`.
#' @param type LAMMPS type.
#' @return A data.frame with columns `atom`, `ql`, `qlbar`.
#' @export
steinhardt <- function(frame, order_l = 6L, cutoff = 3.5, type = 0L,
                      graph = "mutual") {
  UseMethod("steinhardt")
}

#' @export
steinhardt.dseams_frame <- function(frame, order_l = 6L, cutoff = 3.5, type = 0L,
                                    graph = "mutual") {
  out <- steinhardt_(frame$ptr, cutoff, as.integer(type), as.integer(order_l), graph)
  class(out) <- c("dseams_steinhardt", "data.frame")
  out
}

#' Local q4 and q6 on the cutoff shell.
#'
#' Both degrees use the solid-harmonic kernel. `q4bar` and `q6bar` are the
#' Lechner-Dellago averages.
#'
#' @param frame A `dseams_frame`.
#' @param cutoff Neighbour cutoff in angstrom.
#' @param type LAMMPS type. 0 uses the first atom's type.
#' @return A data.frame with columns `atom`, `q4`, `q6`, `q4bar`, `q6bar`.
#' @export
q4_q6 <- function(frame, cutoff = 3.5, type = 0L) {
  UseMethod("q4_q6")
}

#' @export
q4_q6.dseams_frame <- function(frame, cutoff = 3.5, type = 0L) {
  out <- q4_q6_(frame$ptr, cutoff, as.integer(type))
  class(out) <- c("dseams_q4_q6", "data.frame")
  out
}
