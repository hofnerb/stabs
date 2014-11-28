\name{stabsel.stabsel}
\alias{stabsel.stabsel}
\title{
    Change Parameters of Stability Selection
}
\description{
  Method to change the parameters \code{cutoff}, \code{PFER} and
  \code{assumption} of stability selection that can be altered without
  the need to re-run the subsampling process.
}
\usage{
\method{stabsel}{stabsel}(x, cutoff, PFER, assumption = x$assumption, ...)
}
\arguments{
  \item{x}{an object that results from a call to \code{\link{stabsel}}.}
  \item{cutoff}{cutoff between 0.5 and 1. Preferably a value between 0.6
    and 0.9 should be used.}
  \item{PFER}{upper bound for the per-family error rate. This
    specifies the amount of falsely selected base-learners, which is
    tolerated. See details.}
  \item{assumption}{ Defines the type of assumptions on the
    distributions of the selection probabilities and simultaneous
    selection probabilities. Only applicable for
    \code{sampling.type = "SS"}. For \code{sampling.type = "MB"} we
    always use code{"none"}.}
  \item{\dots}{ additional arguments that are currently ignored. }
}
\details{

  This function allows to alter the parameters \code{cutoff},
  \code{PFER} and \code{assumption} of a fitted stability selection
  result. All other parameters are re-used from the original stability
  selection results. The missing paramter is computed and the selected
  variables are updated accordingly.
}
\value{
  An object of class \code{\link{stabsel}}. For details see there.
}
\seealso{
  \code{\link{stabsel}} for the generic function,
  \code{\link{stabsel_parameters}} for the computation of error bounds,
  \code{\link{fitfun}} for available fitting functions and
  \code{\link{plot.stabsel}} for available plot functions
}
\examples{
  ## make data set available
  data("bodyfat", package = "TH.data")
  ## set seed
  set.seed(1234)

  ####################################################################
  ### using stability selection with Lasso methods:

  if (require("lars")) {
      (stab.lasso <- stabsel(x = bodyfat[, -2], y = bodyfat[,2],
                             fitfun = lars.lasso, cutoff = 0.75,
                             PFER = 1))

      par(mfrow = c(2, 1))
      plot(stab.lasso)

      ## now change the PFER and the assumption:
      (stab.lasso_cf0.93_rconc <- stabsel(stab.lasso, cutoff = 0.93,
                                          assumption = "r-concave"))
      plot(stab.lasso_cf0.93_rconc)
      ## the cutoff did change and hence the PFER and the selected
      ## variables
  }
}
\keyword{nonparametric}
