library("stabs")
library("TH.data")

## load data set
data("bodyfat", package = "TH.data")
## and set seed
set.seed(1234)

if (require("mboost")) {
    ############################################################################
    ### stabsel with mboost
    mod <- glmboost(DEXfat ~ ., data = bodyfat)

    ## compute cutoff ahead of running stabsel to see if it is a sensible
    ## parameter choice.
    ##   p = ncol(bodyfat) - 1 (= Outcome) + 1 ( = Intercept)
    stabsel_parameters(q = 3, PFER = 1, p = ncol(bodyfat) - 1 + 1,
                       sampling.type = "MB")
    ## the same:
    stabsel_parameters(mod, q = 3, PFER = 1, sampling.type = "MB")
    ## the same:
    stabsel(mod, q = 3, PFER = 1, sampling.type = "MB", eval = FALSE)

    ## now run stability selection
    (sbody <- stabsel(mod, q = 3, PFER = 1, sampling.type = "MB"))
    opar <- par(mai = par("mai") * c(1, 1, 1, 2.7))
    plot(sbody, type = "paths", ymargin = 6)
    par(opar)
    plot(sbody)
}

if (require("glmnet")) {
    ############################################################################
    ### run stability selection with lasso (from glmnet)
    stab <- stabsel(x = as.matrix(bodyfat[, -2]), y = bodyfat[,2],
                    cutoff = 0.75, PFER = 1, fitfun = glmnet.lasso)
    stab
    plot(stab, type = "maxsel")
    plot(stab, type = "paths")

    ### compare results with hdi
    if (require("hdi")) {
        stab_hdi <- stability(x = as.matrix(bodyfat[, -2]), y = bodyfat[,2],
                              EV = 1, threshold = 0.75)
        stab_hdi
        sort(stab_hdi$freq)
    }

    stab <- stabsel(x = as.matrix(bodyfat[, -2]), y = bodyfat[,2],
                    cutoff = 0.75, PFER = 1, sampling.type = "MB",
                    fitfun = glmnet.lasso)
    stab
    plot(stab, type = "path")
    plot(stab, type = "maxsel")
}

if (require("lars")) {
    ############################################################################
    ### with package lars
    set.seed(1234)
    stab.stepwise <- stabsel(x = as.matrix(bodyfat[, -2]), y = bodyfat[,2],
                             fitfun = lars.stepwise,
                             cutoff = 0.75, PFER = 1, sampling.type = "MB")
    stab.stepwise
    plot(stab.stepwise, type = "maxsel")

    set.seed(1234)
    stab <- stabsel(x = as.matrix(bodyfat[, -2]), y = bodyfat[,2],
                    fitfun = lars.lasso,
                    cutoff = 0.75, PFER = 1, sampling.type = "MB")
    stab
    plot(stab, type = "maxsel")

    ## check data.frame interface
    set.seed(1234)
    stab.df <- stabsel(x = bodyfat[, -2], y = bodyfat[,2],
                       fitfun = lars.lasso,
                       cutoff = 0.75, PFER = 1, sampling.type = "MB")
    stab.df
    stopifnot(all.equal(stab$max, stab.df$max))
    ## and with explicit intercept?
    set.seed(1234)
    stab.int <- stabsel(x = bodyfat[, -2], y = bodyfat[,2], intercept = TRUE,
                        fitfun = lars.lasso,
                        cutoff = 0.75, PFER = 1, sampling.type = "MB")
    stab.int
    stopifnot(all.equal(stab$max, stab.int$max[-1]))

    ############################################################################
    ### use args.fitfun
    set.seed(1234)
    stab.args <- stabsel(x = bodyfat[, -2], y = bodyfat[,2],
                         fitfun = lars.lasso,
                         args.fitfun = list(type = "stepwise"),
                         cutoff = 0.75, PFER = 1, sampling.type = "MB")
    stab.args
    stopifnot(all.equal(stab.stepwise$max, stab.args$max))
}

if (require("mboost")) {
    ############################################################################
    ### get length of formula
    fm <- DEXfat ~ bbs(age) + bbs(waistcirc) + bols(hipcirc) + bols(elbowbreath)
    length(strsplit(deparse(fm), " \\+ ")[[1]])

    ## now this depends on the data...
    fm <- DEXfat ~ .
    length(strsplit(deparse(fm), " \\+ ")[[1]])
}


################################################################################
### check if phat and max are OK

if (require("glmnet")) {
    ### lasso (from glmnet)
    stab <- stabsel(x = as.matrix(bodyfat[, -2]), y = bodyfat[,2],
                    cutoff = 0.75, PFER = 1, fitfun = glmnet.lasso)
    stopifnot(all.equal(stab$max, stab$phat[, ncol(stab$phat)]))
}

if (require("mboost")) {
    ### mboost
    stab <- stabsel(mod, q = 3, PFER = 1, sampling.type = "MB")
    stopifnot(all.equal(stab$max, stab$phat[, ncol(stab$phat)]))
}

if (require("lars")) {
    ### lasso (from lars)
    stab <- stabsel(x = as.matrix(bodyfat[, -2]), y = bodyfat[,2],
                    fitfun = lars.lasso,
                    cutoff = 0.75, PFER = 1, sampling.type = "MB")
    stopifnot(all.equal(stab$max, stab$phat[, ncol(stab$phat)]))

    ### stepwise (from lars)
    stab <- stabsel(x = as.matrix(bodyfat[, -2]), y = bodyfat[,2],
                    fitfun = lars.stepwise,
                    cutoff = 0.75, PFER = 1, sampling.type = "MB")
    stopifnot(all.equal(stab$max, stab$phat[, ncol(stab$phat)]))

    ## what if phat is not available?
    lars.lasso2 <- function(x, y, q, ...) {
        if (!require("lars"))
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
        ret <- logical(ncol(x))
        ret[selected] <- TRUE
        names(ret) <- colnames(x)
        ## return both
        return(list(selected = ret))
    }
    stab <- stabsel(x = as.matrix(bodyfat[, -2]), y = bodyfat[,2],
                    fitfun = lars.lasso2,
                    cutoff = 0.75, PFER = 1, sampling.type = "MB")
    plot(stab, type = "paths")
    ## works.
}


### check if stabsel continues even if a fold leads to an error
if (require("lars")) {
    ############################################################################
    ### with package lars

    check_folds <- function(folds, B, n, sampling.type) {
        if (sampling.type == "SS") {
            folds <- cbind(folds, rep(1, n) - folds)
        }
        folds
    }

    ## set modified check_folds function as function in stabs
    env <- as.environment("package:stabs")
    env2 <- getNamespace("stabs")
    ## store original function
    check_folds_orig <- env[["check_folds"]]
    ## reset it in namespace and package
    unlockBinding("check_folds", env)
    unlockBinding("check_folds", env2)
    env[["check_folds"]] <- check_folds
    env2[["check_folds"]] <- check_folds

    set.seed(1234)
    ## define folds that produce an error
    folds_cc <- folds <- subsample(rep(1, nrow(bodyfat)), B = 50)
    folds[1, 1] <- NA

    stab1 <- stabsel(x = as.matrix(bodyfat[, -2]), y = bodyfat[,2],
                     fitfun = lars.lasso,
                     cutoff = 0.75, PFER = 1, sampling.type = "SS",
                     folds = folds)

    stab2 <- stabsel(x = as.matrix(bodyfat[, -2]), y = bodyfat[,2],
                     fitfun = lars.lasso,
                     cutoff = 0.75, PFER = 1, sampling.type = "SS",
                     folds = folds, papply = lapply)

    stab <- stabsel(x = as.matrix(bodyfat[, -2]), y = bodyfat[,2],
                    fitfun = lars.lasso,
                    cutoff = 0.75, PFER = 1, sampling.type = "SS",
                    folds = folds_cc, papply = lapply)
    stab1
    stab2
    stab

    env[["check_folds"]] <- check_folds_orig
    env2[["check_folds"]] <- check_folds_orig
}

if (require("mboost")) {
    ############################################################################
    ### with package lars

    check_folds <- function(folds, B, n, sampling.type) {
        if (sampling.type == "SS") {
            folds <- cbind(folds, rep(1, n) - folds)
        }
        folds
    }

    ## set modified check_folds function as function in stabs
    env <- as.environment("package:stabs")
    env2 <- getNamespace("stabs")
    ## store original function
    check_folds_orig <- env[["check_folds"]]
    ## reset it in namespace and package
    unlockBinding("check_folds", env)
    unlockBinding("check_folds", env2)
    env[["check_folds"]] <- check_folds
    env2[["check_folds"]] <- check_folds

    set.seed(1234)
    ## define folds that produce an error
    folds_cc <- folds <- subsample(rep(1, nrow(bodyfat)), B = 50)
    folds[1, 1] <- NA

    mod <- glmboost(DEXfat ~ ., data = bodyfat)
    stab1 <- stabsel(mod, cutoff = 0.75, PFER = 1, sampling.type = "SS",
                     folds = folds)

    stab2 <- stabsel(mod, cutoff = 0.75, PFER = 1, sampling.type = "SS",
                     folds = folds, papply = lapply)

    stab <- stabsel(mod, cutoff = 0.75, PFER = 1, sampling.type = "SS",
                    folds = folds_cc, papply = lapply)
    stab1
    stab2
    stab

    env[["check_folds"]] <- check_folds_orig
    env2[["check_folds"]] <- check_folds_orig
}
