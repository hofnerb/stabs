library("stabs")
library("mboost")
library("glmnet")
library("lars")
library("TH.data")

## load data set
data("bodyfat", package = "TH.data")
## and set seed
set.seed(1234)

################################################################################
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

################################################################################
### run stability selection with lasso (from glmnet)
stab <- stabsel(x = as.matrix(bodyfat[, -2]), y = bodyfat[,2],
                cutoff = 0.75, PFER = 1)
stab
plot(stab, type = "maxsel")

### compare results with hdi
if (require("hdi")) {
    stab_hdi <- stability(x = as.matrix(bodyfat[, -2]), y = bodyfat[,2],
                          EV = 1, threshold = 0.75)
    stab_hdi
    sort(stab_hdi$freq)
}

stab <- stabsel(x = as.matrix(bodyfat[, -2]), y = bodyfat[,2],
                cutoff = 0.75, PFER = 1, sampling.type = "MB")
stab
plot(stab, type = "maxsel")

################################################################################
### with package lars
stab <- stabsel(x = as.matrix(bodyfat[, -2]), y = bodyfat[,2],
                fitfun = lars.stepwise,
                cutoff = 0.75, PFER = 1, sampling.type = "MB")
stab
plot(stab, type = "maxsel")

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

################################################################################
if (FALSE) {
## Lasso
library("lars")

lars <- lars(x = as.matrix(bodyfat[, -2]), y = bodyfat[,2], max.steps = 2)
coef(lars)[10,]
lars <- update(lars, x = as.matrix(bodyfat[1:10, -2]), y = bodyfat[1:10,2])

library("glmnet")
mod <- glmnet(x = as.matrix(bodyfat[, -2]), y = bodyfat[,2],
              lambda = seq(0, 20, length = 1000))
plot(mod)
cf <- coef(mod)[, apply(coef(mod), 2, function(x) sum(x != 0)) == 3]
selected <- names(cf[cf[, 1] != 0, 1])

## caret: http://topepo.github.io/caret/modelList.html
}
