\name{parameters}
\alias{parameters}
\alias{stabsel_parameters.stabsel}

\title{ Method to Extract Parameters}
\description{
  Extract stability selection parameters, i.e., tuning parameters, from
  a \code{stabsel} object.
}
\usage{
parameters(object)

## extract parameters from a stabsel model
## (same as parameters(p) )
\method{stabsel_parameters}{stabsel}(p, ...)
}
\arguments{

  \item{object}{ an object of class \code{"stabsel"} or
    \code{"stabsel_parameters"}. }
  \item{p}{an object of class \code{"stabsel"}.}
  \item{...}{ additional arguments, currently not used.}
}
\value{
  An object of class \code{\link{stabsel_parameters}} with a special
  \code{print} method. See there for details.
}
\seealso{
  \code{\link{stabsel}} to run stability selection and
  \code{\link{stabsel_parameters}} for details on the parameters.
}
\keyword{methods}
