\name{subsample}
\alias{subsample}
\title{ Draw Random Subsamples }
\description{
  Set up weight matrix for subsampling with sample proportion \eqn{1/2}
  to be used with \code{\link{stabsel}}.
}
\usage{
subsample(weights, B = 100, strata = NULL)
}
\arguments{
  \item{weights}{ a numeric vector of weights for the model to be cross-validated.}
  \item{B}{ number of folds, per default 25 for \code{bootstrap} and
    \code{subsampling} and 10 for \code{kfold}.}
  \item{strata}{ a factor of the same length as \code{weights} for
    stratification.}
}
\details{
  The function \code{subsample} can be used to build an appropriate
  weight matrix to be used with \code{\link{stabsel}}. See there for
  more details.

  If \code{strata} is defined sampling is performed in each stratum
  separately thus preserving the distribution of the \code{strata}
  variable in each fold.
}
\seealso{
  \code{\link{stabsel}}
}
\examples{
  ## just a low-dimensional example
  subsample(weights = rep(1, 10), B = 50)
}
\keyword{nonparametric}

