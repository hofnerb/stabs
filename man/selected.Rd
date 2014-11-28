\name{selected}
\alias{selected}
\alias{selected.stabsel}

\title{ Method to Extract Selected Variables }
\description{
  Extract selected variables from a \code{stabsel} object.
}
\usage{
selected(object, ...)
\method{selected}{stabsel}(object, ...)
}
\arguments{
  \item{object}{ an object of class \code{"stabsel"}. }
  \item{\dots}{ additional arguments passed to specific \code{selected}
    methods. }
}
\details{
  The ids of variables selected during the stability selection process
  can be extracted using \code{selected()}.
}
\keyword{methods}
