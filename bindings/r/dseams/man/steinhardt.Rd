\name{steinhardt}
\alias{steinhardt}
\alias{steinhardt.dseams_frame}
\title{Steinhardt order parameters}
\usage{
steinhardt(frame, order_l = 6L, cutoff = 3.5, type = 0L)
}
\arguments{
\item{frame}{A \code{dseams_frame}.}
\item{order_l}{Spherical-harmonic degree. 3, 4, 6, 8, and 12 are computed.
Degrees 4 and 6 on an orthorhombic or restricted-triclinic cell use the
solid-harmonic kernel.}
\item{cutoff}{Neighbour cutoff in angstrom. An orthorhombic or
restricted-triclinic cell uses the cutoff shell. Any other cell keeps the
mutual four nearest inside \code{cutoff + 1.5}.}
\item{type}{LAMMPS type.}
}
\value{
A data.frame with columns \code{atom}, \code{ql}, \code{qlbar}.
}
\description{
Local Steinhardt \eqn{q_l} and the Lechner-Dellago average \eqn{\bar q_l}.
}
