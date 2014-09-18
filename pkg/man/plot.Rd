\name{plot.stabsel}
\alias{plot}
\alias{plot.stabsel}
\title{
    Plotting Method for Stability Selection
}
\description{
  Display results of stability selection.
}
\usage{
\method{plot}{stabsel}(x, main = deparse(x$call), type = c("maxsel", "paths"),
     col = NULL, ymargin = 10, np = sum(x$max > 0),
     labels = NULL, ...)
}
\arguments{
  \item{x}{object of class \code{stabsel}.}
  \item{main}{main title for the plot.}
  \item{type}{plot type; either stability paths (\code{"paths"}) or a
    plot of the maximum selection frequency (\code{"maxsel"}).}
  \item{col}{a vector of colors; Typically, one can specify a single
     color or one color for each variable. Per default, colors depend on
     the maximal selection frequency of the variable and range from grey
     to red.}
  \item{ymargin}{(temporarily) specifies the y margin of of the plot in
    lines (see argument \code{"mar"} of function \code{\link{par}}).
    This only affects the right margin for \code{type = "paths"} and
    the left margin for \code{type = "maxsel"}. Explicit user specified
    margins are kept and are not overwritten.}
  \item{np}{number of variables to plot for the maximum selection
    frequency plot (\code{type = "maxsel"}); the first \code{np}
    variables with highest selection frequency are plotted.}
  \item{labels}{variable labels for the plot; one label per base-learner
    must be specified. Per default, names of base-learners are used.}
  \item{\dots}{additional arguments to plot functions.}
}
\details{

  This function implements the stability selection procedure
  by Meinshausen and Buehlmann (2010) and the improved error bounds by
  Shah and Samworth (2013).

  Two of the three arguments \code{cutoff}, \code{q} and \code{PFER}
  \emph{must} be specified. The per-family error rate (PFER), i.e., the
  expected number of false positives \eqn{E(V)}, where \eqn{V} is the
  number of false positives, is bounded by the argument \code{PFER}.

  As controlling the PFER is more conservative as controlling the
  family-wise error rate (FWER), the procedure also controlls the FWER,
  i.e., the probability of selecting at least one non-influential
  variable (or model component) is less than \code{PFER}.

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
  \item{sampling.type}{the sampling type used for stability selection.}
  \item{assumption}{the assumptions made on the selection
    probabilities.}
  \item{call}{the call.}
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
  \code{\link{stabsel}}
}
\examples{

  library("mboost")
  data("bodyfat", package = "TH.data")

  ### low-dimensional example
  mod <- glmboost(DEXfat ~ ., data = bodyfat)

  ## compute cutoff ahead of running stabsel to see if it is a sensible
  ## parameter choice.
  ##   p = ncol(bodyfat) - 1 (= Outcome) + 1 ( = Intercept)
  stabsel_parameters(q = 3, PFER = 1, p = ncol(bodyfat) - 1 + 1,
                     sampling.type = "MB")
  ## the same:
  stabsel(mod, q = 3, PFER = 1, sampling.type = "MB", eval = FALSE)

  ## now run stability selection; to make results reproducible
  set.seed(1234)
  (sbody <- stabsel(mod, q = 3, PFER = 1, sampling.type = "MB"))
  opar <- par(mai = par("mai") * c(1, 1, 1, 2.7))
  plot(sbody)
  par(opar)

  plot(sbody, type = "maxsel", ymargin = 6)
}
\keyword{nonparametric}
