## Regression test for stability selection with QUIC
context("QUIC stability selection regression test")
library(stabs)
library(huge)
set.seed(10010)
dat.hubs <- huge.generator(n=N, d=40, graph="hub")
pcutoff <- 0.75
PFER <- 10
s.hubs <- stabsel(x=dat.hubs$data, fitfun=stabs.quic, 
                  cutoff=pcutoff, PFER=PFER, verbose=FALSE)

load("quic.hub.test.Rda")
testthat::test_that("Selected same", expect_equal(s.hubs$selected, s.hubs.orig$selected))
