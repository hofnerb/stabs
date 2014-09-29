## Generic implementation of stability selection
stabsel <- function(x, ...) {
    UseMethod("stabsel", x)
}


### TODO: parallelization ala cvrisk needed! Seems to work?
### TODO: Use same arguments for .mboost und .formula
### TODO: Should y be a matrix? Perhaps we need this for survival data which
### might be specified as a matrix?
stabsel.matrix <- function(x, y, fitfun = lars.lasso, args.fitfun = list(),
                           cutoff, q, PFER,
                           folds = subsample(rep(1, nrow(x)), B = B),
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
    fit_model <- function(i, folds, q, args.fitfun) {
        inbag <- as.logical(folds[, i])
        do.call(fitfun, c(list(x = x[inbag, ], y = y[inbag, ], q = q),
                          args.fitfun))
    }

    nms <- colnames(x)
    ret <- run_stabsel(fitter = fit_model, args.fitter = args.fitfun,
                n = n, p = p, cutoff = cutoff, q = q,
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


### TODO: What about weights?
### n <- sum(weights)
##    stabsel.formula <- function(formula, data, weights = rep(1, nrow(data)), fitfun = glmnet.lasso,
##                                args.fitfun = list(), p = NULL, cutoff, q, PFER,
##                                folds = subsample(rep(1, nrow(x)), B = B),
##                                B = ifelse(sampling.type == "MB", 100, 50),
##                                assumption = c("unimodal", "r-concave", "none"),
##                                sampling.type = c("SS", "MB"),
##                                papply = mclapply, verbose = TRUE, FWER, eval = TRUE,
##                                ...) {
##
##        warning("This function is very experimental at the moment")
##
##        cll <- match.call()
##        ## TODO: ??? How do I get this for all formulae?
##        ## perhaps one can fit the model once and obtain p <- length(coef)
##
##        ## try to guess p
##        if (is.null(p))
##            p <- length(strsplit(deparse(formula), " \\+ ")[[1]])
##        n <- nrow(data)
##
##        ## needed here to make B and folds happy
##        sampling.type <- match.arg(sampling.type)
##        if (sampling.type == "MB")
##            assumption <- "none"
##        else
##            assumption <- match.arg(assumption)
##
##        ## define fitting function;
##        ## the function implicitly knows formula and data as it is defined in this environment
##        fit_model <- function(i, folds, q, args.fitfun) {
##            inbag <- as.logical(folds[, i])
##            do.call(fitfun, c(list(formula = formula, data = data, q = q),
##                              args.fitfun))
##        }
##
##        nms <- colnames(x)
##        ret <- run_stabsel(fitter = fit_model, args.fitter = args.fitfun,
##                           n = n, p = p, cutoff = cutoff, q = q,
##                           PFER = PFER, folds = folds, B = B, assumption = assumption,
##                           sampling.type = sampling.type, papply = papply,
##                           verbose = verbose, FWER = FWER, eval = eval, names = nms, ...)
##        ret$call <- cll
##        ret$call[[1]] <- as.name("stabsel")
##        return(ret)
##
##    }

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

