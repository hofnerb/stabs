\name{stabsel}
\alias{stabsel}
\alias{stabsel.matrix}
\alias{stabsel.data.frame}
\title{
    Stability Selection
}
\description{
    Selection of influential variables or model components with error control.
}
\usage{
## generic stability selection funcion
stabsel(x, ...)

## a method to fit models with stability selection
\method{stabsel}{matrix}(x, y, fitfun = lars.lasso,
        args.fitfun = list(), cutoff, q, PFER,
        folds = subsample(rep(1, nrow(x)), B = B),
        B = ifelse(sampling.type == "MB", 100, 50),
        assumption = c("unimodal", "r-concave", "none"),
        sampling.type = c("SS", "MB"),
        papply = mclapply, verbose = TRUE, FWER, eval = TRUE, ...)

## essentially a wrapper for data.frames (see details)
\method{stabsel}{data.frame}(x,  y, intercept = FALSE, ...)
}
\arguments{
  \item{x}{a \code{\link{matrix}} or a \code{\link{data.frame}}
    containing the predictors.}
  \item{y}{a vector or matrix containing the outcome.}
  \item{intercept}{logical. If \code{x} is a \code{\link{data.frame}},
    this argument determines if the resulting model matrix should
    contain a separate intercept or not.}
  \item{fitfun}{a function that takes the arguments \code{x}, \code{y}
    as above, and additionally the number of variables to include in
    each model \code{q}. The function then needs to fit the model and to
    return a logical vector that indicates which variable was selected
    (among the \code{q} selected variables).}
  \item{args.fitfun}{a named list containing additional arguments that are
    passed to the fitting function; see also argument \code{args} in
    \code{\link{do.call}}.}
  \item{cutoff}{cutoff between 0.5 and 1. Preferably a value between 0.6
    and 0.9 should be used.}
  \item{q}{number of (unique) selected variables (or groups of variables
    depending on the model) that are selected on each subsample.}
  \item{PFER}{upper bound for the per-family error rate. This
    specifies the amount of falsely selected base-learners, which is
    tolerated. See details.}
  \item{folds}{ a weight matrix with number of rows equal to the number
    of observations, see \code{\link{subsample}}. Usually one should not
    change the default here as subsampling with a fraction of \eqn{1/2}
    is needed for the error bounds to hold. One usage scenario where
    specifying the folds by hand might be the case when one has
    dependent data (e.g. clusters) and thus wants to draw clusters
    (i.e., multiple rows together) not individuals.}
  \item{assumption}{ Defines the type of assumptions on the
    distributions of the selection probabilities and simultaneous
    selection probabilities. Only applicable for
    \code{sampling.type = "SS"}. For \code{sampling.type = "MB"} we
    always use code{"none"}.}
  \item{sampling.type}{ use sampling scheme of of Shah & Samworth
    (2013), i.e., with complementarty pairs (\code{sampling.type = "SS"}),
    or the original sampling scheme of Meinshausen & Buehlmann (2010).}
  \item{B}{ number of subsampling replicates. Per default, we use 50
    complementary pairs for the error bounds of Shah & Samworth (2013)
    and 100 for the error bound derived in  Meinshausen & Buehlmann
    (2010). As we use \eqn{B} complementray pairs in the former case
    this leads to \eqn{2B} subsamples.}
  \item{papply}{ (parallel) apply function, defaults to
    \code{\link[parallel]{mclapply}}. Alternatively, \code{parLapply}
    can be used. In the latter case, usually more setup is needed (see
    example of \code{\link[mboost]{cvrisk}} for some details).}
  \item{verbose}{ logical (default: \code{TRUE}) that determines wether
    \code{warnings} should be issued. }
  \item{FWER}{ deprecated. Only for compatibility with older versions,
    use PFER instead.}
  \item{eval}{ logical. Determines whether stability selection is
    evaluated (\code{eval = TRUE}; default) or if only the parameter
    combination is returned.}
  \item{\dots}{ additional arguments to parallel apply methods such as
    \code{\link{mclapply}}.}
}
\details{

  This function implements the stability selection procedure by
  Meinshausen and Buehlmann (2010) and the improved error bounds by Shah
  and Samworth (2013). For details see also Hofner et al. (2014). The
  error bounds are implemented in the function
  \code{\link{stabsel_parameters}}. Two of the three arguments
  \code{cutoff}, \code{q} and \code{PFER} \emph{must} be specified. The
  per-family error rate (PFER), i.e., the expected number of false
  positives \eqn{E(V)}, where \eqn{V} is the number of false positives,
  is bounded by the argument \code{PFER}.

  As controlling the PFER is more conservative as controlling the
  family-wise error rate (FWER), the procedure also controlls the FWER,
  i.e., the probability of selecting at least one non-influential
  variable (or model component) is less than \code{PFER}.

  Predefined \code{\link{fitfuns}} functions exist but more can be
  easily implemented. Note that stepwise regression methods are usually
  not advised as they tend to be relatively unstable. See example below.

  The function \code{stabsel} for \code{\link{data.frame}}s is
  essentially just a wrapper to the \code{\link{matrix}} function with
  the same argments. The only difference is that in a pre-processing
  step, the data set is converted to a model matrix using the function
  \code{\link{model.matrix}}. The additional argument \code{intercept}
  determines if an explicit intercept should be added to the model
  matrix. This is often not neccessary but depends on the \code{fitfun}.

}
\value{
  An object of class \code{stabsel} with a special \code{print} method.
  The object has the following elements:
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
  \item{call}{the call.}
}
\references{

  B. Hofner, L. Boccuto and M. Goeker (2014), Controlling false
  discoveries in high-dimensional situations: Boosting with stability
  selection. \emph{Technical Report}, arXiv:1411.1285.\cr
  \url{http://arxiv.org/abs/1411.1285}.

  N. Meinshausen and P. Buehlmann (2010), Stability selection.
  \emph{Journal of the Royal Statistical Society, Series B},
  \bold{72}, 417--473.

  R.D. Shah and R.J. Samworth (2013), Variable selection with error
  control: another look at stability selection. \emph{Journal of the Royal
  Statistical Society, Series B}, \bold{75}, 55--80.

}
\seealso{
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
      (stab.stepwise <- stabsel(x = bodyfat[, -2], y = bodyfat[,2],
                                fitfun = lars.stepwise, cutoff = 0.75,
                                PFER = 1))
      par(mfrow = c(2, 1))
      plot(stab.lasso, main = "Lasso")
      plot(stab.stepwise, main = "Stepwise Selection")
      ## --> stepwise selection seems to be quite unstable even in this low
      ##     dimensional example!
  }

  ####################################################################
  ### using stability selection directly on computed boosting models
  ### from mboost

  if (require("mboost")) {
      ### low-dimensional example
      mod <- glmboost(DEXfat ~ ., data = bodyfat)

      ## compute cutoff ahead of running stabsel to see if it is a sensible
      ## parameter choice.
      ##   p = ncol(bodyfat) - 1 (= Outcome) + 1 ( = Intercept)
      stabsel_parameters(q = 3, PFER = 1, p = ncol(bodyfat) - 1 + 1,
                         sampling.type = "MB")
      ## the same:
      stabsel(mod, q = 3, PFER = 1, sampling.type = "MB", eval = FALSE)

      ## now run stability selection
      (sbody <- stabsel(mod, q = 3, PFER = 1, sampling.type = "MB"))
      opar <- par(mai = par("mai") * c(1, 1, 1, 2.7))
      plot(sbody)
      par(opar)

      plot(sbody, type = "maxsel", ymargin = 6)
  }
}
\keyword{nonparametric}
