################################################################################
## fit functions
##
## functions need to take arguments x, y and q and return a logical vector that
## indicates which variable was selected
##
################################################################################

glmnet.lasso <- function(x, y, q, type = c("conservative", "anticonservative"), ...) {
    if (!requireNamespace("glmnet", quietly=TRUE))
        stop("Package ", sQuote("glmnet"), " needed but not available")
    
    if (is.data.frame(x)) {
        message("Note: ", sQuote("x"),
                " is coerced to a model matrix without intercept")
        x <- model.matrix(~ . - 1, x)
    }
    
    if ("lambda" %in% names(list(...)))
        stop("It is not permitted to specify the penalty parameter ", sQuote("lambda"),
             " for lasso when used with stability selection.")
    
    ## fit model
    type <- match.arg(type)
    if (type == "conservative")
        fit <- glmnet::glmnet(x, y, pmax = q, ...)
    if (type == "anticonservative")
        fit <- glmnet::glmnet(x, y, dfmax = q - 1, ...)
    
    ## which coefficients are non-zero?
    selected <- predict(fit, type = "nonzero")
    selected <- selected[[length(selected)]]
    ret <- logical(ncol(x))
    ret[selected] <- TRUE
    names(ret) <- colnames(x)
    ## compute selection paths
    cf <- fit$beta
    sequence <- as.matrix(cf != 0)
    ## return both
    return(list(selected = ret, path = sequence))
}

lars.lasso <- function(x, y, q, ...) {
    if (!requireNamespace("lars", quietly=TRUE))
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

lars.stepwise <- function(x, y, q, ...) {
    if (!requireNamespace("lars", quietly=TRUE))
        stop("Package ", sQuote("lars"), " needed but not available")

    if (is.data.frame(x)) {
        message("Note: ", sQuote("x"),
                " is coerced to a model matrix without intercept")
        x <- model.matrix(~ . - 1, x)
    }

    ## fit model
    fit <- lars::lars(x, y, max.steps = q, type = "stepwise", ...)

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

glmnet.lasso_maxCoef <- function(x, y, q, ...) {
    if (!requireNamespace("glmnet", quietly=TRUE))
        stop("Package ", sQuote("glmnet"), " needed but not available")

    if (is.data.frame(x)) {
        message("Note: ", sQuote("x"),
                " is coerced to a model matrix without intercept")
        x <- model.matrix(~ . - 1, x)
    }

    args <- list(...)
    if (!("lambda" %in% names(args)) && length(args$lambda) != 1)
        stop("Please specify a fixed (!) value of ", sQuote("lambda"),
             ", which is small enough that at least ", sQuote("q"),
             " variables can be selected.")

    ## fit model
    fit <- glmnet::glmnet(x, y, ...)

    ## which coefficients are the q biggest
    selected <- order(coef(fit)[-1])[1:q]
    ret <- logical(ncol(x))
    ret[selected] <- TRUE
    names(ret) <- colnames(x)
    ## return selection
    return(list(selected = ret, path = NULL))
}


## mboost.glmboost <- function(formula, data, weights, ...) {
##     if (!require("mboost"))
##         stop("Package ", sQuote("mboost"), " needed but not available")
##
##     ## fit model
##     fit <- glmboost(formula, data, weights, ...)
## }
