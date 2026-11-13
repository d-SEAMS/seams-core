\name{primitive_rings}
\alias{primitive_rings}
\alias{primitive_rings.dseams_frame}
\title{Primitive ring census}
\usage{
primitive_rings(frame, cutoff = 3.5, type = 0L, graph = "mutual",
                max_depth = 6L)
}
\arguments{
\item{frame}{A \code{dseams_frame}.}
\item{cutoff}{Neighbour cutoff in angstrom.}
\item{type}{LAMMPS type.}
\item{graph}{\code{"mutual"} or \code{"cutoff"}.}
\item{max_depth}{Largest ring retained.}
}
\value{
A data.frame with columns \code{size} and \code{n}.
}
\description{
Franzblau shortest-path rings on the bond graph.
On Ice Ih with \code{graph = "cutoff"} the six-ring row is 256.
See \code{\link{chill_plus}} for the example that asserts it.
}
