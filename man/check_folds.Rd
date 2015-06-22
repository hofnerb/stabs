\name{check_folds}
\alias{check_folds}
\title{
  Check if folds result from subsampling with p = 0.5.
}
\description{
  (Internal) function that checks if folds result from subsampling with
  p = 0.5 and adds complementary pairs if needed.
}
\usage{
check_folds(folds, B, n, sampling.type)
}
\arguments{
  \item{folds}{ a weight matrix that represents the subsamples.}
  \item{B}{ number of subsampling replicates. }
  \item{n}{ the number of observations; needed for internal checks.}
  \item{sampling.type}{ sampling type to be used.}
}
\details{

  This is an internal function used to check if \code{folds} are
  specified correctly. For details (e.g. on arguments) see
  \code{\link{stabsel}}.
}
\value{
  A matrix containing the \code{folds}, possibly after adding the
  complementary pairs.
}
\references{

  B. Hofner, L. Boccuto and M. Goeker (2015), Controlling false
  discoveries in high-dimensional situations: Boosting with stability
  selection. \emph{BMC Bioinformatics}, 16:144.\cr
  doi:10.1186/s12859-015-0575-3.

}
\seealso{
  For details see \code{\link{stabsel}}.
}
\keyword{helper}
