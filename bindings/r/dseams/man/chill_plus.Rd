\name{chill_plus}
\alias{chill_plus}
\alias{chill_plus.dseams_frame}
\title{CHILL+ ice labels}
\usage{
chill_plus(frame, cutoff = 3.5, type = 0L, graph = "cutoff",
          rule = "chill+")
}
\arguments{
\item{frame}{A \code{dseams_frame}.}
\item{cutoff}{Neighbour cutoff in angstrom. On the mutual graph this is
the cell-size hint.}
\item{type}{LAMMPS type. 0 uses the first atom.}
\item{graph}{\code{"cutoff"} (the default, every pair closer than
\code{cutoff}) or \code{"mutual"} (the mutual four nearest neighbours
inside \code{cutoff + 1.5}).}
\item{rule}{\code{"chill+"} is the published rule: exactly four neighbours
inside the cutoff, q3 over the four nearest, and the Table 1 neighbour
test. \code{"ovito-chill+"} averages q3 over the whole cutoff shell and
skips that neighbour test.}
}
\value{
A data.frame with columns \code{atom}, \code{type}, \code{label},
\code{name}. The \code{counts} attribute is the nine-label histogram.
}
\description{
Per-atom CHILL+ labels from libyodaLib.
Ice Ih (\code{input/traj/ice1h.lammpstrj}, every atom) is 128 hexagonal
waters. The example asserts that count and the 256 primitive six-rings.
}
\examples{
root <- Sys.getenv("DSEAMS_ROOT", "")
if (nzchar(root)) {
  frame <- read_frame(file.path(root, "input", "traj", "ice1h.lammpstrj"),
                      type = -1L)
  ice <- chill_plus(frame, graph = "cutoff")
  stopifnot(sum(ice$name == "hexagonal") == 128L)
  rings <- primitive_rings(frame, graph = "cutoff", max_depth = 6L)
  stopifnot(rings$n[rings$size == 6L] == 256L)
}
}
