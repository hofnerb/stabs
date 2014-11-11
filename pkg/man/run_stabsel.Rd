\name{run_stabsel}
\alias{run_stabsel}
\title{
  Run Stability Selection
}
\description{
  (Internal) function that is used to run stability selection (i.e. to
  apply the fit-function to the subsamples. This function is not intended
  to be directly called.
}
\usage{
run_stabsel(fitter, args.fitter, n, p, cutoff, q, PFER, folds, B, assumption,
            sampling.type, papply, verbose, FWER, eval, names, ...)
}
\arguments{
  \item{fitter}{a function to fit the model on subsamples. See argument
    \code{fitfun} of \code{\link{stabsel}} for details.}
  \item{args.fitter}{a named list containing additional arguments that are
    passed to \code{fitter}. See argument \code{args.fitfun}
    \code{\link{stabsel}} for details.}
  \item{n}{ the number of observations; needed for internal checks.}
  \item{p}{ number of possible predictors (including intercept if
    applicable).}
  \item{cutoff}{cutoff between 0.5 and 1. }
  \item{q}{number of (unique) selected variables (or groups of variables
    depending on the model) that are selected on each subsample.}
  \item{PFER}{ upper bound for the per-family error rate.}
  \item{folds}{ a weight matrix that represents the subsamples.}
  \item{B}{ number of subsampling replicates. }
  \item{assumption}{ distributional assumption.}
  \item{sampling.type}{ sampling type to be used.}
  \item{papply}{ (parallel) apply function.}
  \item{verbose}{ logical (default: \code{TRUE}) that determines wether
    \code{warnings} should be issued. }
  \item{FWER}{ deprecated. Only for compatibility with older versions,
    use PFER instead.}
  \item{eval}{ logical. Determines whether stability selection is
    evaluated.}
  \item{names}{ variable names that are used to label the results.}
  \item{\dots}{ additional arguments to be passed to next function.}
}
\details{

  This is an internal function that fits the actual models to the
  subsamples, i.e., this is the work horse that runs stability
  selection. Usually, one should use  \code{\link{stabsel}}, which
  internally calls \code{run_stabsel}.

  \code{run_stabsel} can be used by expert users to implement stability
  selection methods for new model types.

  For details (e.g. on arguments) see \code{\link{stabsel}}.
}
\value{
  An object of class \code{stabsel} with the following elements:
   \item{phat}{selection probabilities.}
  \item{selected}{elements with maximal selection probability greater
    \code{cutoff}.}
  \item{max}{maximum of selection probabilities.}
  \item{cutoff}{cutoff used.}
  \item{q}{average number of selected variables used.}
  \item{PFER}{per-family error rate.}
  \item{p}{the number of effects subject to selection.}
  \item{sampling.type}{the sampling type used for stability selection.}
  \item{assumption}{the assumptions made on the selection
    probabilities.}
}
\references{

  B. Hofner, L. Boccuto and M. Goeker (2014), Controlling false
  discoveries in high-dimensional situations: Boosting with stability
  selection. \emph{Technical Report}, arXiv:1411.1285.\cr
  \url{http://arxiv.org/abs/1411.1285}.

}
\seealso{
  For details see \code{\link{stabsel}}.
}
\keyword{helper}
