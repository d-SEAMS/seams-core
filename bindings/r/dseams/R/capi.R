# Copyright (c) d-SEAMS developers
# .Call wrappers over include/seams.h. The libseams_lang frame API stays in dseams.R.

#' @export
seams_system_new <- function() {
  .Call(`_dseams_seams_r_system`)
}

#' @export
seams_system_read <- function(sys, path, frame = 1L, type = -1L) {
  .Call(`_dseams_seams_r_read`, sys, normalizePath(path, mustWork = TRUE), as.integer(frame), as.integer(type))
}

#' @export
seams_system_count <- function(sys) {
  .Call(`_dseams_seams_r_count`, sys)
}

#' @export
seams_chill <- function(sys, cutoff = 3.5, type = 0L, graph = 0L, rule = 0L) {
  .Call(`_dseams_seams_r_chill`, sys, cutoff, as.integer(type), as.integer(graph), as.integer(rule))
}

#' @export
seams_steinhardt <- function(sys, order, cutoff = 3.5, average = FALSE) {
  .Call(`_dseams_seams_r_steinhardt`, sys, cutoff, as.integer(order), as.integer(average))
}

#' @export
seams_cages <- function(sys, cutoff = 3.5, graph = 2L, k = 4L, complete = FALSE) {
  .Call(`_dseams_seams_r_cages`, sys, cutoff, as.integer(graph), as.integer(k), as.integer(complete))
}

#' @export
seams_signature_cages <- function(sys, spec, cutoff = 3.5, graph = 0L) {
  .Call(`_dseams_seams_r_signature`, sys, cutoff, as.integer(graph), spec)
}

#' @export
seams_neighbours_knn <- function(sys, k = 4L, cutoff = 3.5, mutual = TRUE) {
  .Call(`_dseams_seams_r_knn`, sys, as.integer(k), cutoff, as.integer(mutual))
}

#' @export
seams_rdf <- function(sys, rmax = 6, bins = 8L) {
  .Call(`_dseams_seams_r_rdf`, sys, rmax, as.integer(bins))
}

#' @export
seams_running_cn <- function(sys, rmax = 6, bins = 8L) {
  .Call(`_dseams_seams_r_running_cn`, sys, rmax, as.integer(bins))
}

#' @export
seams_coordination <- function(sys, rmax = 3.5, bins = 8L) {
  .Call(`_dseams_seams_r_coordination`, sys, rmax, as.integer(bins))
}

#' @export
seams_hbonds <- function(sys, cutoff = 3.5, oxygen_type = 1L, hydrogen_type = 2L,
                         dist = 2.42, angle = 30) {
  .Call(`_dseams_seams_r_hbonds`, sys, cutoff, as.integer(oxygen_type),
        as.integer(hydrogen_type), dist, angle)
}

#' @export
seams_density_axis <- function(sys, type = 0L, bins = 4L, axis = 2L) {
  .Call(`_dseams_seams_r_density`, sys, as.integer(type), as.integer(bins), as.integer(axis))
}

#' @export
seams_fingerprint <- function(sys, cutoff = 3.5, hops = 2L) {
  .Call(`_dseams_seams_r_fingerprint`, sys, cutoff, as.integer(hops))
}

#' @export
seams_contact_pairs <- function(sys, spec) {
  .Call(`_dseams_seams_r_pairs`, sys, spec)
}

#' @export
seams_domains <- function(sys, spec, cutoff = 3.5, subset = 0L) {
  .Call(`_dseams_seams_r_domains`, sys, cutoff, spec, as.integer(subset))
}

#' @export
seams_prisms <- function(sys, cutoff = 3.5) {
  .Call(`_dseams_seams_r_prisms`, sys, cutoff)
}

#' @export
seams_monolayer <- function(sys, cutoff = 3.5, max_depth = 6L, sheet_area = 100) {
  .Call(`_dseams_seams_r_monolayer`, sys, cutoff, as.integer(max_depth), sheet_area)
}

#' @export
seams_rings_vertices <- function(sys, cutoff = 3.5, max_depth = 6L, index = 0L) {
  .Call(`_dseams_seams_r_vertices`, sys, cutoff, as.integer(max_depth), as.integer(index))
}

#' @export
seams_site_table <- function(sys, spec) {
  .Call(`_dseams_seams_r_site`, sys, spec)
}
