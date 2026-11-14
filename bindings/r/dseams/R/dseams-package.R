#' @keywords internal
#' @useDynLib dseams, .registration = TRUE
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  n <- getOption("dseams.threads")
  if (is.null(n)) {
    limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(limit) && !identical(toupper(limit), "FALSE")) {
      n <- min(2L, get_threads_())
    }
  }
  if (!is.null(n)) {
    set_threads_(as.integer(n))
  }
}
