\name{jump_tau90}
\alias{jump_tau90}
\title{Jump-rotor tau_90}
\usage{
jump_tau90(first, second, dt = 1, oxygen_type = 1L, hydrogen_type = 2L)
}
\arguments{
\item{first}{The earlier \code{dseams_frame}.}
\item{second}{The later \code{dseams_frame}.}
\item{dt}{Time between the frames.}
\item{oxygen_type}{LAMMPS type of oxygen.}
\item{hydrogen_type}{LAMMPS type of hydrogen.}
}
\value{
\code{dt} when any water's H-H vector rotates by at least 90 degrees,
and \code{NaN} otherwise.
}
\description{
Jump-rotor tau_90 between two frames.
}
