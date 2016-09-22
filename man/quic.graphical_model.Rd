\name{quic.graphical_model}
\alias{quic.graphical_model}
\alias{getLamPath}

\title{Stability selection fit function for sparse inverse covariance using QUIC}
\usage{
quic.graphical_model(x, y, q, ...)
getLamPath(max, min, len, log = FALSE)
}
\arguments{
\item{x}{data matrix}
\item{y}{data matrix}
\item{q}{number of variables}
\item{...}{extra arguments, currently only \code{lams} (optional) a lambda path}
\item{max}{maximum value for regularization (lambda)}
\item{min}{min value for lambda}
\item{len}{length of path}
\item{log}{log spacing}
}
\value{
list with selected variables and path.
}
\description{
Stability selection fit function for sparse inverse covariance using QUIC.
}
\author{
  Code contributed by Richard Beare.
}
\details{
This is a wrapper for QUIC to be used in stability selection. Pass it as
the fit function to \code{\link{stabsel}}. If you supply your own graphical model fitter, please
set the class to \code{"graphical_model"}.

\code{getLamPath} can be used to create a regularization path.
}
\examples{
if (require(huge) & require(QUIC)) {
 set.seed(10010)
 dat.hubs <- huge.generator(n=1000, d=40, graph="hub")
 pcutoff <- 0.75
 PFER <- 10
 s.hubs <- stabsel(x=dat.hubs$data, fitfun=quic.graphical_model, 
                  cutoff=pcutoff, PFER=PFER)

}
}
\seealso{
\code{\link{stabsel}}
}

