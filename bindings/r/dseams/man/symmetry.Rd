\name{symmetry}
\alias{symmetry}
\alias{symmetry.dseams_frame}
\title{Shell symmetry label}
\usage{
symmetry(frame, shell_rmsd = 0.12)
}
\arguments{
\item{frame}{A \code{dseams_frame}.}
\item{shell_rmsd}{Shell RMSD tolerance.}
}
\value{
Integer codes: 0 ice Ih, 1 ice Ic, 2 stacking-disordered, -1 otherwise.
}
\description{
Per-atom shell symmetry label.
}
