################################################################################
## fit functions
##
## functions need to take arguments x, y and q and return a logical vector that
## indicates which variable was selected
##
################################################################################

glmnet.lasso <- function(x, y, q, ...) {
    if (!require("glmnet"))
        stop("Package ", sQuote("glmnet"), " needed but not available")

    if (is.data.frame(x)) {
        message("Note: ", sQuote("x"),
                " is coerced to a model matrix without intercept")
        x <- model.matrix(~ . - 1, x)
    }

    ## fit model
    fit <- glmnet(x, y, dfmax = q - 1, ...)

    ## which coefficients are non-zero?
    selected <- predict(fit, type = "nonzero")
    selected <- selected[[length(selected)]]
    ret <- logical(ncol(x))
    ret[selected] <- TRUE
    names(ret) <- colnames(x)
    ret
}

lars.lasso <- function(x, y, q, ...) {
    if (!require("lars"))
        stop("Package ", sQuote("lars"), " needed but not available")

    if (is.data.frame(x)) {
        message("Note: ", sQuote("x"),
                " is coerced to a model matrix without intercept")
        x <- model.matrix(~ . - 1, x)
    }

    ## fit model
    fit <- lars(x, y, max.steps = q, ...)

    ## which coefficients are non-zero?
    selected <- unlist(fit$actions)
    ret <- logical(ncol(x))
    ret[selected] <- TRUE
    names(ret) <- colnames(x)
    ret
}

lars.stepwise <- function(x, y, q, ...) {
    if (!require("lars"))
        stop("Package ", sQuote("lars"), " needed but not available")

    if (is.data.frame(x)) {
        message("Note: ", sQuote("x"),
                " is coerced to a model matrix without intercept")
        x <- model.matrix(~ . - 1, x)
    }

    ## fit model
    fit <- lars(x, y, max.steps = q, type = "stepwise", ...)

    ## which coefficients are non-zero?
    selected <- unlist(fit$actions)
    ret <- logical(ncol(x))
    ret[selected] <- TRUE
    names(ret) <- colnames(x)
    ret
}

## mboost.glmboost <- function(formula, data, weights, ...) {
##     if (!require("mboost"))
##         stop("Package ", sQuote("mboost"), " needed but not available")
##
##     ## fit model
##     fit <- glmboost(formula, data, weights, ...)
## }
