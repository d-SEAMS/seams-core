\name{chill_plus}
\alias{chill_plus}
\alias{chill_plus.dseams_frame}
\title{CHILL+ ice labels}
\usage{
chill_plus(frame, cutoff = 3.5, type = 0L, graph = "cutoff")
}
\arguments{
\item{frame}{A \code{dseams_frame}.}
\item{cutoff}{Neighbour cutoff in angstrom. On the mutual graph this is
the cell-size hint.}
\item{type}{LAMMPS type. 0 uses the first atom.}
\item{graph}{\code{"cutoff"} (the default, every pair closer than
\code{cutoff}) or \code{"mutual"} (the mutual four nearest neighbours
inside \code{cutoff + 1.5}).}
}
\value{
A data.frame with columns \code{atom}, \code{type}, \code{label},
\code{name}. The \code{counts} attribute is the nine-label histogram.
}
\description{
Per-atom CHILL+ labels from libyodaLib.
}
