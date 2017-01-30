\name{Fitting Functions}
\alias{fitfun}
\alias{fitfuns}
\alias{glmnet.lasso}
\alias{glmnet.lasso_maxCoef}
\alias{lars.lasso}
\alias{lars.stepwise}

\title{
  Fit Functions for Stability Selection
}
\description{
  Functions that fit a model until \eqn{q} variables are selected and
  that returns the indices (and names) of the selected variables.
}
\usage{
## package lars:
lars.lasso(x, y, q, ...)
lars.stepwise(x, y, q, ...)

## package glmnet:
glmnet.lasso(x, y, q, type = c("conservative", "anticonservative"), ...)
glmnet.lasso_maxCoef(x, y, q, ...)
}
\arguments{
  \item{x}{a matrix containing the predictors or an object of class
    \code{"mboost"}.}
  \item{y}{a vector or matrix containing the outcome.}
  \item{q}{number of (unique) selected variables (or groups of variables
    depending on the model) that are selected on each subsample.}
  \item{type}{ a charachter vector specifying if the number of selected
    variables per subsample is \eqn{\leq q}{<= q} (\code{type = "conservative"}) 
    or \eqn{\geq q}{>= q} (\code{type = "anticonservative"}). The conservative
    version \emph{ensures} that the PFER is controlled.}
  \item{\dots}{ additional arguments passed to the underlying fitting function.
    See the example on \code{glmnet.lasso_maxCoef} in \code{\link{stabsel}} for 
    the specification of additional arguments via \code{\link{stabsel}}.}
}
\details{
  All fitting functions are named after the package and the type of
  model that is fitted: \code{package_name.model}, e.g.,
  \code{glmnet.lasso} stands for a lasso model that is fitted using the
  package \pkg{glmnet}.

  \code{glmnet.lasso_maxCoef} fits a lasso model with a given penalty parameter 
  and returns the \code{q} largest coefficients. If one wants to use
  \code{glmnet.lasso_maxCoef}, one must specify the penalty parameter
  \code{lambda} (via the \code{\dots} argument) or in
  \code{\link{stabsel}} via \code{args.fitfun(lambda = )}. Note that usually, 
  the penalty parameter cannot be specified but is chosen such that \code{q}
  variables are selected. For an example on how to use 
  \code{glmnet.lasso_maxCoef} see \code{\link{stabsel}}.
  
}
\value{
  A named list with elements
  \item{selected}{logical. A vector that indicates which variable was
    selected.}
  \item{path}{logical. A matrix that indicates which variable was
    selected in which step. Each row represents one variable, the
    columns represent the steps.}
}
\seealso{
\code{\link{stabsel}} for stability selection itself, and \code{\link{quic.graphical_model}} for stability selection for graphical models.
}
\examples{
  if (require("TH.data")) {
      ## make data set available
      data("bodyfat", package = "TH.data")
  } else {
      ## simulate some data if TH.data not available. 
      ## Note that results are non-sense with this data.
      bodyfat <- matrix(rnorm(720), nrow = 72, ncol = 10)
  }
  
  if (require("lars")) {
      ## selected variables
      lars.lasso(bodyfat[, -2], bodyfat[,2], q = 3)$selected
      lars.stepwise(bodyfat[, -2], bodyfat[,2], q = 3)$selected
  }
  
  if (require("glmnet")) {
      glmnet.lasso(bodyfat[, -2], bodyfat[,2], q = 3)$selected
      ## selection path
      glmnet.lasso(bodyfat[, -2], bodyfat[,2], q = 3)$path
  
      ## Using the anticonservative glmnet.lasso (see args.fitfun):
      stab.glmnet <- stabsel(x = bodyfat[, -2], y = bodyfat[,2],
                             fitfun = glmnet.lasso, 
                             args.fitfun = list(type = "anticonservative"), 
                             cutoff = 0.75, PFER = 1)
  }
}
\keyword{models}
\keyword{nonlinear}
\keyword{nonparametric}
