## Generic implementation of stability selection
stabsel <- function(x, ...) {
    UseMethod("stabsel", x)
}


### TODO: parallelization ala cvrisk needed! Seems to work?
### TODO: Use same arguments for .mboost und .formula
### TODO: Add args.fitfun = list() that takes additional arguments for the
### fitter.
### TODO: Should y be a matrix? Perhaps we need this for survival data which
### might be specified as a matrix?
stabsel.matrix <- function(x, y, fitfun = glmnet.lasso, cutoff, q, PFER,
                           folds = cv(rep(1, nrow(x)), type = "subsampling", B = B),
                           B = ifelse(sampling.type == "MB", 100, 50),
                           assumption = c("unimodal", "r-concave", "none"),
                           sampling.type = c("SS", "MB"),
                           papply = mclapply, verbose = TRUE, FWER, eval = TRUE,
                           ...) {

    cll <- match.call()
    p <- ncol(x) ## TODO: what about intercept?
    n <- nrow(x)

    ## needed here to make B and folds happy
    sampling.type <- match.arg(sampling.type)
    if (sampling.type == "MB")
        assumption <- "none"
    else
        assumption <- match.arg(assumption)

    ## make sure that y is a matrix
    if (!is.matrix(y))
        y <- matrix(y, ncol = 1)

    if (n != nrow(y))
        stop(sQuote("x"), " and ", sQuote("y"),
             " must have the same number of observations")

    ## define fitting function;
    ## the function implicitly knows x and y as it is defined in this environment
    fit_model <- function(i, folds, q) {
        inbag <- as.logical(folds[, i])
        do.call(fitfun, list(x = x[inbag, ], y = y[inbag, ], q = q))
    }

    nms <- colnames(x)
    ret <- run_stabsel(fitter = fit_model, n = n, p = p, cutoff = cutoff, q = q,
                PFER = PFER, folds = folds, B = B, assumption = assumption,
                sampling.type = sampling.type, papply = papply,
                verbose = verbose, FWER = FWER, eval = eval, names = nms, ...)
    ret$call <- cll
    ret$call[[1]] <- as.name("stabsel")
    return(ret)
}

stabsel.data.frame <- function(x, y, intercept = FALSE, ...) {
    if (intercept) {
        x <- model.matrix(~ ., x)
    } else {
        x <- model.matrix(~ . - 1, x)
    }
    stabsel(x, y, ...)
}

stabsel.formula <- function(formula, ...) {
}

stabsel.mboost <- function(x, cutoff, q, PFER,
                    folds = cv(model.weights(x), type = "subsampling", B = B),
                    B = ifelse(sampling.type == "MB", 100, 50),
                    assumption = c("unimodal", "r-concave", "none"),
                    sampling.type = c("SS", "MB"),
                    papply = mclapply, verbose = TRUE, FWER, eval = TRUE, ...) {

    cll <- match.call()
    p <- length(variable.names(x))
    ibase <- 1:p

    sampling.type <- match.arg(sampling.type)
    if (sampling.type == "MB")
        assumption <- "none"
    else
        assumption <- match.arg(assumption)

    B <- ncol(folds)

    pars <- stabsel_parameters(p = p, cutoff = cutoff, q = q,
                               PFER = PFER, B = B,
                               verbose = verbose, sampling.type = sampling.type,
                               assumption = assumption, FWER = FWER)
    ## return parameter combination only if eval == FALSE
    if (!eval)
        return(pars)

    cutoff <- pars$cutoff
    q <- pars$q
    PFER <- pars$PFER

    fun <- function(model) {
        xs <- selected(model)
        qq <- sapply(1:length(xs), function(x) length(unique(xs[1:x])))
        xs[qq > q] <- xs[1]
        xs
    }
    if (sampling.type == "SS") {
        ## use complementary pairs
        folds <- cbind(folds, model.weights(x) - folds)
    }
    ss <- cvrisk(x, fun = fun,
                 folds = folds,
                 papply = papply, ...)

    if (verbose){
        qq <- sapply(ss, function(x) length(unique(x)))
        sum_of_violations <- sum(qq < q)
        if (sum_of_violations > 0)
            warning(sQuote("mstop"), " too small in ",
                    sum_of_violations, " of the ", ncol(folds),
                    " subsampling replicates to select ", sQuote("q"),
                    " base-learners; Increase ", sQuote("mstop"),
                    " bevor applying ", sQuote("stabsel"))
    }


    ## if grid specified in '...'
    if (length(list(...)) >= 1 && "grid" %in% names(list(...))) {
        m <- max(list(...)$grid)
    } else {
        m <- mstop(x)
    }
    ret <- matrix(0, nrow = length(ibase), ncol = m)
    for (i in 1:length(ss)) {
        tmp <- sapply(ibase, function(x)
            ifelse(x %in% ss[[i]], which(ss[[i]] == x)[1], m + 1))
        ret <- ret + t(sapply(tmp, function(x) c(rep(0, x - 1), rep(1, m - x + 1))))
    }

    phat <- ret / length(ss)
    rownames(phat) <- names(variable.names(x))
    if (extends(class(x), "glmboost"))
        rownames(phat) <- variable.names(x)
    ret <- list(phat = phat, selected = which((mm <- apply(phat, 1, max)) >= cutoff),
                max = mm, cutoff = cutoff, q = q, PFER = PFER,
                sampling.type = sampling.type, assumption = assumption,
                call = cll)
    ret$call[[1]] <- as.name("stabsel")
    class(ret) <- "stabsel"
    ret
}

## stabsel.default <- function(x, cutoff, q, PFER,
##                     folds = cv(model.weights(object), type = "subsampling",
##                                B = ifelse(sampling.type == "MB", 100, 50)),
##                     assumption = c("unimodal", "r-concave", "none"),
##                     sampling.type = c("SS", "MB"),
##                     papply = mclapply, verbose = TRUE, FWER, eval = TRUE, ...) {
##
##     call <- match.call()
##     p <- length(variable.names(object))
##     ibase <- 1:p
##
##     sampling.type <- match.arg(sampling.type)
##     if (sampling.type == "MB")
##         assumption <- "none"
##     else
##         assumption <- match.arg(assumption)
##
##     B <- ncol(folds)
##
##     pars <- stabsel_parameters(p = p, cutoff = cutoff, q = q,
##                                PFER = PFER, B = B,
##                                verbose = verbose, sampling.type = sampling.type,
##                                assumption = assumption, FWER = FWER)
##
##     ## return parameter combination only if eval == FALSE
##     if (!eval)
##         return(pars)
##
##     cutoff <- pars$cutoff
##     q <- pars$q
##     PFER <- pars$PFER
##
##     fun <- function(model) {
##         xs <- selected(model)
##         qq <- sapply(1:length(xs), function(x) length(unique(xs[1:x])))
##         xs[qq > q] <- xs[1]
##         xs
##     }
##     if (sampling.type == "SS") {
##         ## use complementary pairs
##         folds <- cbind(folds, model.weights(object) - folds)
##     }
##     ss <- cvrisk(object, fun = fun,
##                  folds = folds,
##                  papply = papply, ...)
##
##     if (verbose){
##         qq <- sapply(ss, function(x) length(unique(x)))
##         sum_of_violations <- sum(qq < q)
##         if (sum_of_violations > 0)
##             warning(sQuote("mstop"), " too small in ",
##                     sum_of_violations, " of the ", ncol(folds),
##                     " subsampling replicates to select ", sQuote("q"),
##                     " base-learners; Increase ", sQuote("mstop"),
##                     " bevor applying ", sQuote("stabsel"))
##     }
##
##
##     ## if grid specified in '...'
##     if (length(list(...)) >= 1 && "grid" %in% names(list(...))) {
##         m <- max(list(...)$grid)
##     } else {
##         m <- mstop(object)
##     }
##     ret <- matrix(0, nrow = length(ibase), ncol = m)
##     for (i in 1:length(ss)) {
##         tmp <- sapply(ibase, function(x)
##             ifelse(x %in% ss[[i]], which(ss[[i]] == x)[1], m + 1))
##         ret <- ret + t(sapply(tmp, function(x) c(rep(0, x - 1), rep(1, m - x + 1))))
##     }
##
##     phat <- ret / length(ss)
##     rownames(phat) <- names(variable.names(object))
##     if (extends(class(object), "glmboost"))
##         rownames(phat) <- variable.names(object)
##     ret <- list(phat = phat, selected = which((mm <- apply(phat, 1, max)) >= cutoff),
##                 max = mm, cutoff = cutoff, q = q, PFER = PFER,
##                 sampling.type = sampling.type, assumption = assumption,
##                 call = call)
##     class(ret) <- "stabsel"
##     ret
## }

stabsel_parameters <- function(p, ...)
    UseMethod("stabsel_parameters")

stabsel_parameters.default <- function(p, cutoff, q, PFER,
                               B = ifelse(sampling.type == "MB", 100, 50),
                               assumption = c("unimodal", "r-concave", "none"),
                               sampling.type = c("SS", "MB"),
                               verbose = FALSE, FWER, ...) {

    sampling.type <- match.arg(sampling.type)
    if (sampling.type == "MB")
        assumption <- "none"
    else
        assumption <- match.arg(assumption)


    ## only two of the four arguments can be specified
    if ((nmiss <- sum(missing(PFER), missing(cutoff),
                      missing(q), missing(FWER))) != 2) {
        if (nmiss > 2)
            stop("Two of the three argumnets ",
                 sQuote("PFER"), ", ", sQuote("cutoff"), " and ", sQuote("q"),
                 " must be specifed")
        if (nmiss < 2)
            stop("Only two of the three argumnets ",
                 sQuote("PFER"), ", ", sQuote("cutoff"), " and ", sQuote("q"),
                 " can be specifed at the same time")
    }

    if (!missing(FWER)) {
        if (!missing(PFER))
            stop(sQuote("FWER"), " and ", sQuote("PFER"),
                 " cannot be spefified at the same time")
        PFER <- FWER
        warning(sQuote("FWER"), " is deprecated. Use ", sQuote("PFER"),
                " instead.")
    }

    if ((!missing(PFER) || !missing(FWER)) && PFER < 0)
        stop(sQuote("PFER"), " must be greater 0")

    if (!missing(cutoff) && (cutoff < 0.5 | cutoff > 1))
        stop(sQuote("cutoff"), " must be between 0.5 and 1")

    if (!missing(q)) {
        if (p < q)
            stop("Average number of selected base-learners ", sQuote("q"),
                 " must be smaller \n  than the number of base-learners",
                 " specified in the model ", sQuote("object"))
        if (q < 0)
            stop("Average number of selected base-learners ", sQuote("q"),
                 " must be greater 0")
    }

    if (missing(cutoff)) {
        if (assumption == "none") {
            cutoff <- min(1, tmp <- (q^2 / (PFER * p) + 1) / 2)
            upperbound <- q^2 / p / (2 * cutoff - 1)
        } else {
            if (assumption == "unimodal") {
                cutoff <- tmp <- optimal_cutoff(p, q, PFER, B,
                                                assumption = assumption)
                upperbound <- q^2 / p / um_const(cutoff, B, theta = q/p)
            } else {
                cutoff <- tmp <- optimal_cutoff(p, q, PFER, B,
                                                assumption = assumption)
                upperbound <- minD(q, p, cutoff, B) * p
            }
        }
        upperbound <- signif(upperbound, 3)
        if (verbose && tmp > 0.9 && upperbound - PFER > PFER/2) {
            warning("Upper bound for PFER > ", PFER,
                    " for the given value of ", sQuote("q"),
                    " (true upper bound = ", round(upperbound, 2), ")")
        }
    }

    if (missing(q)) {
        if (assumption == "none") {
            q <- floor(sqrt(PFER * (2 * cutoff - 1) * p))
            upperbound <- q^2 / p / (2 * cutoff - 1)
        } else {
            if (assumption == "unimodal") {
                q <- optimal_q(p, cutoff, PFER, B, assumption = assumption)
                upperbound <- q^2 / p / um_const(cutoff, B, theta = q/p)
            } else {
                q <- optimal_q(p, cutoff, PFER, B, assumption = assumption)
                upperbound <- minD(q, p, cutoff, B) * p
            }
        }
        upperbound <- signif(upperbound, 3)
        if (verbose && upperbound - PFER > PFER/2)
            warning("Upper bound for PFER > ", PFER,
                    " for the given value of ", sQuote("cutoff"),
                    " (true upper bound = ", upperbound, ")")
    }

    if (missing(PFER)) {
        if (assumption == "none") {
            upperbound <- PFER <- q^2 / p / (2 * cutoff - 1)
        } else {
            if (assumption == "unimodal") {
                upperbound <- PFER <- q^2 / p / um_const(cutoff, B, theta = q/p)
            } else {
                upperbound <- PFER <- minD(q, p, cutoff, B) * p
            }
        }
        upperbound <- signif(upperbound, 3)
    }

    if (verbose && PFER >= p)
        warning("Upper bound for PFER larger than the number of base-learners.")

    res <- list(cutoff = cutoff, q = q, PFER = upperbound,
                sampling.type = sampling.type, assumption = assumption)
    class(res) <- "stabsel_parameters"
    res
}

stabsel_parameters.mboost <- function(p, ...) {
    stabsel(p, ..., eval = FALSE)
}
