\name{stabsel_parameters}
\alias{stabsel_parameters}
\alias{stabsel_parameters.default}
\title{
  Compute Error Bounds for Stability Selection
}
\description{
  Compute
}
\usage{
## generic function to compute missing parameter from the other two parameters;
## (internally used within stabsel)
stabsel_parameters(p, ...)

\method{stabsel_parameters}{default}(p, cutoff, q, PFER,
                   B = ifelse(sampling.type == "MB", 100, 50),
                   assumption = c("unimodal", "r-concave", "none"),
                   sampling.type = c("SS", "MB"),
                   verbose = FALSE, FWER, ...)
}
\arguments{
  \item{p}{ number of possible predictors (including intercept if
    applicable).}
  \item{cutoff}{cutoff between 0.5 and 1. Preferably a value between 0.6
    and 0.9 should be used.}
  \item{q}{number of (unique) selected variables (or groups of variables
    depending on the model) that are selected on each subsample.}
  \item{PFER}{upper bound for the per-family error rate. This
    specifies the amount of falsely selected base-learners, which is
    tolerated. See details.}
  \item{B}{ number of subsampling replicates. Per default, we use 50
    complementary pairs for the error bounds of Shah & Samworth (2013)
    and 100 for the error bound derived in  Meinshausen & Buehlmann
    (2010). As we use \eqn{B} complementray pairs in the former case
    this leads to \eqn{2B} subsamples.}
  \item{assumption}{ Defines the type of assumptions on the
    distributions of the selection probabilities and simultaneous
    selection probabilities. Only applicable for
    \code{sampling.type = "SS"}. For \code{sampling.type = "MB"} we
    always use code{"none"}.}
  \item{sampling.type}{ use sampling scheme of of Shah & Samworth
    (2013), i.e., with complementarty pairs (\code{sampling.type = "SS"}),
    or the original sampling scheme of Meinshausen & Buehlmann (2010).}
  \item{verbose}{ logical (default: \code{TRUE}) that determines wether
    \code{warnings} should be issued. }
  \item{FWER}{ deprecated. Only for compatibility with older versions,
    use PFER instead.}
  \item{\dots}{ additional arguments to be passed to next function.}
}
\details{

  This function implements the error bounds for stability selection
  by Meinshausen and Buehlmann (2010) and the improved error bounds by
  Shah and Samworth (2013).

  Two of the three arguments \code{cutoff}, \code{q} and \code{PFER}
  \emph{must} be specified. The per-family error rate (PFER), i.e., the
  expected number of false positives \eqn{E(V)}, where \eqn{V} is the
  number of false positives, is bounded by the argument \code{PFER}.

  For more details see also \code{\link{stabsel}}.
}
\value{
  An object of class \code{stabsel_parameters} with a special \code{print} method.
  The object has the following elements:
  \item{cutoff}{cutoff used.}
  \item{q}{average number of selected variables used.}
  \item{PFER}{per-family error rate.}
  \item{sampling.type}{the sampling type used for stability selection.}
  \item{assumption}{the assumptions made on the selection
    probabilities.}
}
\references{

  N. Meinshausen and P. Buehlmann (2010), Stability selection.
  \emph{Journal of the Royal Statistical Society, Series B},
  \bold{72}, 417--473.

  R.D. Shah and R.J. Samworth (2013), Variable selection with error
  control: another look at stability selection. \emph{Journal of the Royal
  Statistical Society, Series B}, \bold{75}, 55--80.

}
\seealso{
  For more details see also \code{\link{stabsel}}.
}
\keyword{helper}
