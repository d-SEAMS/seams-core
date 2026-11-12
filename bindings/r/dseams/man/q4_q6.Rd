\name{q4_q6}
\alias{q4_q6}
\alias{q4_q6.dseams_frame}
\title{Local q4 and q6}
\usage{
q4_q6(frame, cutoff = 3.5, type = 0L)
}
\arguments{
\item{frame}{A \code{dseams_frame}.}

\item{cutoff}{Neighbour cutoff in angstrom.}

\item{type}{LAMMPS type. 0 uses the first atom's type.}
}
\value{
A data.frame with columns \code{atom}, \code{q4}, \code{q6}, \code{q4bar},
\code{q6bar}.
}
\description{
Local Steinhardt \eqn{q_4} and \eqn{q_6} on the cutoff shell, from the
solid-harmonic kernel, and the Lechner-Dellago averages.
}
