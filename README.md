stabs
=====

[![Build Status](https://travis-ci.org/hofnerb/stabs.svg)](https://travis-ci.org/hofnerb/stabs)
[![Build status](https://ci.appveyor.com/api/projects/status/tlo7dbrevje1f2du?svg=true)](https://ci.appveyor.com/project/hofnerb/stabs)
[![Coverage Status](https://coveralls.io/repos/hofnerb/stabs/badge.svg?branch=master&service=github)](https://coveralls.io/github/hofnerb/stabs?branch=master)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/stabs)](https://cran.r-project.org/package=stabs)
[![](http://cranlogs.r-pkg.org/badges/stabs)](https://cran.r-project.org/package=stabs)

`stabs` implements resampling procedures to assess the stability of selected
variables with additional finite sample error control for high-dimensional
variable selection procedures such as Lasso or boosting. Both, standard
stability selection (Meinshausen & Bühlmann, 2010) and complementarty pairs
stability selection with improved error bounds (Shah & Samworth, 2013) are
implemented. The package can be combined with arbitrary user specified variable
selection approaches.

For an expanded and executable version of this file please see
```r
vignette("Using_stabs", package = "stabs")
```

## Installation

- Current version (from CRAN):

        install.packages("stabs")

- Latest development version from GitHub:

        library("devtools")
        install_github("hofnerb/stabs")

  To be able to use the `install_github()` command, one needs to install `devtools` first:

        install.packages("devtools")

## Using stabs

A simple example of how to use `stabs` with package `lars`:

```
library("stabs")
library("lars")
## make data set available
data("bodyfat", package = "TH.data")
## set seed
set.seed(1234)

## lasso
(stab.lasso <- stabsel(x = bodyfat[, -2], y = bodyfat[,2],
                       fitfun = lars.lasso, cutoff = 0.75,
                       PFER = 1))

## stepwise selection
(stab.stepwise <- stabsel(x = bodyfat[, -2], y = bodyfat[,2],
                          fitfun = lars.stepwise, cutoff = 0.75,
                          PFER = 1))

## plot results
par(mfrow = c(2, 1))
plot(stab.lasso, main = "Lasso")
plot(stab.stepwise, main = "Stepwise Selection")
```

We can see that stepwise selection seems to be quite unstable even in this low
dimensional example!

### User-specified variable selection approaches

To use `stabs` with user specified functions, one can specify an own `fitfun`.
These need to take arguments `x` (the predictors), `y` (the outcome) and `q` the
number of selected variables as defined for stability selection. Additional
arguments to the variable selection method can be handled by `...`. In the
function `stabsel()` these can then be specified as a named list which is given
to `args.fitfun`.

The `fitfun` function then needs to return a named list with two elements
`selected` and `path`:
* `selected` is a vector that indicates which variable was selected.
* `path` is a matrix that indicates which variable was selected in which step.
    Each row represents one variable, the columns represent the steps.
The latter is optional and only needed to draw the complete selection paths.

The following example shows how `lars.lasso` is implemented:
```
lars.lasso <- function(x, y, q, ...) {
    if (!requireNamespace("lars"))
        stop("Package ", sQuote("lars"), " needed but not available")

    if (is.data.frame(x)) {
        message("Note: ", sQuote("x"),
                " is coerced to a model matrix without intercept")
        x <- model.matrix(~ . - 1, x)
    }

    ## fit model
    fit <- lars::lars(x, y, max.steps = q, ...)

    ## which coefficients are non-zero?
    selected <- unlist(fit$actions)
	## check if variables are removed again from the active set
    ## and remove these from selected
    if (any(selected < 0)) {
        idx <- which(selected < 0)
        idx <- c(idx, which(selected %in% abs(selected[idx])))
        selected <- selected[-idx]
    }

    ret <- logical(ncol(x))
    ret[selected] <- TRUE
    names(ret) <- colnames(x)
    ## compute selection paths
    cf <- fit$beta
    sequence <- t(cf != 0)
    ## return both
    return(list(selected = ret, path = sequence))
}
```

To see more examples simply print, e.g., `lars.stepwise`, `glmnet.lasso`, or
`glmnet.lasso_maxCoef`. Please contact me if you need help to integrate your
method of choice.

### Using boosting with stability selection

Instead of specifying a fitting function, one can also use `stabsel` directly on
computed boosting models from
[mboost](https://cran.r-project.org/package=mboost).

```
library("stabs")
library("mboost")
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
plot(sbody, type = "paths")
par(opar)

plot(sbody, type = "maxsel", ymargin = 6)
```

## Citation

To cite the package in publications please use
```
citation("stabs")
```

which will currently give you

```
To cite package 'stabs' in publications use:

  Benjamin Hofner and Torsten Hothorn (2015). stabs: Stability
  Selection with Error Control, R package version R package version
  0.6-0, http://CRAN.R-project.org/package=stabs.

  Benjamin Hofner, Luigi Boccuto and Markus Goeker (2015). Controlling
  false discoveries in high-dimensional situations: Boosting with
  stability selection. BMC Bioinformatics, 16:144.
  doi:10.1186/s12859-015-0575-3

Use ‘toBibtex(citation("stabs"))’ to extract BibTeX references.
```

To obtain BibTeX references use

```
toBibtex(citation("stabs"))
```
