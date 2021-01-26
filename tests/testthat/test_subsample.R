## Regression test for subsampling
context("subsample regression tests")

library("stabs")

set.seed(1907)

test_that("Subsample (without stratification) works", {
    weights <- rep(1, 100)
    ss1 <- subsample(weights, B = 10)
    expect_equal(colSums(ss1), rep(floor(length(weights) * 0.5), 10))
    
    weights <- rep(1, 101)
    ss2 <- subsample(weights, B = 20)
    expect_equal(colSums(ss2), rep(floor(length(weights) * 0.5), 20))
})



test_that("Subsample (with stratification) works", {
    
    weights <- rep(1, 100)
    strata1 <- as.factor(sample(1:4, length(weights), replace = TRUE))
    n_stratum1 <- table(strata1)
    ss1 <- subsample(weights, B = 10, strata = strata1)
    expect_equal(colSums(ss1), 
                 rep(sum(floor(n_stratum1 * 0.5)), 10))
    
    strata2 <- as.factor(sample(1:4, length(weights), replace = TRUE))
    n_stratum2 <- table(strata2)
    ss2 <- subsample(weights, B = 10, strata = strata2)
    expect_equal(colSums(ss2), 
                 rep(sum(floor(n_stratum2 * 0.5)), 10))
})

