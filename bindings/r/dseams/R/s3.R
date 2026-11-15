# Copyright (c) d-SEAMS developers
# print and summary methods. The frame printer stays in dseams.R.

#' Summary of one frame.
#'
#' @param object A `dseams_frame`.
#' @param ... Unused.
#' @return `object`, invisibly.
#' @export
summary.dseams_frame <- function(object, ...) {
  cat("dseams frame", object$frame, "from", object$path, "\n", sep = " ")
  cat("nop", nop_(object$ptr), "type filter", object$type, "\n", sep = " ")
  invisible(object)
}

#' Print CHILL+ labels.
#'
#' @param x A `dseams_ice` data frame.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.dseams_ice <- function(x, ...) {
  cat("CHILL+ labels,", nrow(x), "atoms\n")
  base::print(table(x$name))
  invisible(x)
}

#' Summary of CHILL+ labels.
#'
#' @param object A `dseams_ice` data frame.
#' @param ... Unused.
#' @return `object`, invisibly.
#' @export
summary.dseams_ice <- function(object, ...) {
  base::print(table(object$name))
  invisible(object)
}

#' Print cage counts.
#'
#' @param x A `dseams_cages` data frame.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.dseams_cages <- function(x, ...) {
  cat(
    "cages hexagonal", attr(x, "n_ih"),
    "cubic", attr(x, "n_ic"),
    "both", attr(x, "n_both"),
    "water", attr(x, "n_water"),
    "\n"
  )
  invisible(x)
}

#' Summary of cage names.
#'
#' @param object A `dseams_cages` data frame.
#' @param ... Unused.
#' @return `object`, invisibly.
#' @export
summary.dseams_cages <- function(object, ...) {
  base::print(table(object$name))
  invisible(object)
}

#' Print a ring histogram.
#'
#' @param x A `dseams_rings` data frame.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.dseams_rings <- function(x, ...) {
  base::print(x[, c("size", "n")])
  invisible(x)
}

#' Summary of a ring histogram.
#'
#' @param object A `dseams_rings` data frame.
#' @param ... Unused.
#' @return `object`, invisibly.
#' @export
summary.dseams_rings <- function(object, ...) {
  base::print(object[, c("size", "n")])
  invisible(object)
}

#' Print Steinhardt parameters.
#'
#' @param x A `dseams_steinhardt` data frame.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.dseams_steinhardt <- function(x, ...) {
  cat("Steinhardt parameters,", nrow(x), "atoms, mean ql", mean(x$ql), "\n")
  invisible(x)
}

#' Summary of Steinhardt parameters.
#'
#' @param object A `dseams_steinhardt` data frame.
#' @param ... Unused.
#' @return `object`, invisibly.
#' @export
summary.dseams_steinhardt <- function(object, ...) {
  cat("mean ql", mean(object$ql), "mean qlbar", mean(object$qlbar), "\n")
  invisible(object)
}

#' Print q4 and q6.
#'
#' @param x A `dseams_q4_q6` data frame.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.dseams_q4_q6 <- function(x, ...) {
  cat("q4 q6,", nrow(x), "atoms, mean q6", mean(x$q6), "\n")
  invisible(x)
}

#' Summary of q4 and q6.
#'
#' @param object A `dseams_q4_q6` data frame.
#' @param ... Unused.
#' @return `object`, invisibly.
#' @export
summary.dseams_q4_q6 <- function(object, ...) {
  cat("mean q4", mean(object$q4), "mean q6", mean(object$q6), "\n")
  invisible(object)
}

#' Print a system handle.
#'
#' @param x A `dseams_system`.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.dseams_system <- function(x, ...) {
  n <- tryCatch(seams_system_count(x), error = function(e) NA_integer_)
  cat("dseams system", n, "atoms\n")
  invisible(x)
}

#' Summary of a system handle.
#'
#' @param object A `dseams_system`.
#' @param ... Unused.
#' @return `object`, invisibly.
#' @export
summary.dseams_system <- function(object, ...) {
  cat("external pointer for one libseams frame\n")
  base::print(object)
  invisible(object)
}

#' Print CHILL population counts.
#'
#' @param x A `dseams_chill` vector.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.dseams_chill <- function(x, ...) {
  cat("CHILL counts\n")
  base::print(stats::setNames(as.integer(x), names(x)))
  invisible(x)
}

#' Summary of CHILL population counts.
#'
#' @param object A `dseams_chill` vector.
#' @param ... Unused.
#' @return `object`, invisibly.
#' @export
summary.dseams_chill <- function(object, ...) {
  base::print(object)
  invisible(object)
}
