\name{signature_cages}
\alias{signature_cages}
\alias{signature_cages.dseams_frame}
\title{Closed cages for one ring census}
\usage{
signature_cages(frame, spec = "512", cutoff = 3.5, type = 0L,
                graph = "cutoff")
}
\arguments{
\item{frame}{A \code{dseams_frame}.}
\item{spec}{Ring census such as \code{"512"} or \code{"5:12,6:2"}.}
\item{cutoff}{Neighbour cutoff in angstrom.}
\item{type}{LAMMPS type.}
\item{graph}{\code{"cutoff"} or \code{"mutual"}.}
}
\value{
An integer vector \code{c(n_cages, n_atoms)}.
}
\description{
Closed cages whose faces match \code{spec}.
}
