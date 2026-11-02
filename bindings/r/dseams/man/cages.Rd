\name{cages}
\alias{cages}
\alias{cages.dseams_frame}
\title{Cage membership}
\usage{
cages(frame, cutoff = 3.5, type = 0L, k = 4L, graph = "seeded")
}
\arguments{
\item{frame}{A \code{dseams_frame}.}
\item{cutoff}{Neighbour cutoff in angstrom.}
\item{type}{LAMMPS type.}
\item{k}{Neighbours kept by a k-nearest graph.}
\item{graph}{\code{"seeded"}, \code{"mutual"}, \code{"cutoff"}, or
\code{"knn-union"}.}
}
\value{
A data.frame with columns \code{atom}, \code{cage}, \code{name}.
Attributes \code{n_ih}, \code{n_ic}, \code{n_both}, and \code{n_water}
are disjoint and match \code{seams cages}. Per-atom \code{cage} is
0 water, 1 hexagonal only, 2 double-diamond only, 3 both.
}
\description{
Hexagonal-cage and double-diamond membership.
}
