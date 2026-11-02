\name{read_frame}
\alias{read_frame}
\alias{nframes}
\alias{dseams_frame}
\title{Read one trajectory frame}
\usage{
read_frame(path, frame = 1L, type = 0L)
nframes(path)
}
\arguments{
\item{path}{LAMMPS dump or XYZ file.}
\item{frame}{1-based frame index.}
\item{type}{LAMMPS type filter. 0 guesses oxygen (type 2, else type 1).
A negative value keeps every atom.}
}
\value{
\code{read_frame} returns an object of class \code{dseams_frame}.
\code{nframes} returns the number of LAMMPS timesteps.
}
\description{
Read one frame from a trajectory the d-SEAMS engine can open.
}
