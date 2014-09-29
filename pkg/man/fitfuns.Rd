\name{Fitting Functions}
\alias{fitfun}
\alias{fitfuns}
\alias{glmnet.lasso}
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
glmnet.lasso(x, y, q, ...)
}
\arguments{
  \item{x}{a matrix containing the predictors or an object of class
    \code{"mboost"}.}
  \item{y}{a vector or matrix containing the outcome.}
  \item{q}{number of (unique) selected variables (or groups of variables
    depending on the model) that are selected on each subsample.}
  \item{\dots}{ additional arguments to the underlying fitting function.}
}
\details{
  All fitting functions are named after the package and the type of
  model that is fitted: \code{package_name.model}, e.g.,
  \code{glmnet.lasso} stands for a lasso model that is fitted using the
  package \pkg{glmnet}.
}
\value{
  A named list with elements
  \itemize{
    \item \code{selected} logical. A vector that indicates which variable was
    selected.
    \item \code{path} logical. A matrix that indicates which variabkle was
    selected in which step. Each row represents one variable, the
    columns represent the steps.
  }
}
\examples{
data("bodyfat", package = "TH.data")
## selected variables
lars.lasso(bodyfat[, -2], bodyfat[,2], q = 3)$selected
lars.stepwise(bodyfat[, -2], bodyfat[,2], q = 3)$selected
glmnet.lasso(bodyfat[, -2], bodyfat[,2], q = 3)$selected
}
\keyword{models}
\keyword{nonlinear}
\keyword{nonparametric}
