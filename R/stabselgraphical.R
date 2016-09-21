#' Create a regularization path - copied from pulsar.
#' @export
#' @param max maximum value for regularization (lambda)
#' @param min min value for lambda
#' @param len length of path
#' @param log log spacing
getLamPath <- function (max, min, len, log = FALSE) 
{
  if (max < min) 
    stop("Did you flip min and max?")
  if (log) {
    min <- log(min)
    max <- log(max)
  }
  lams <- seq(max, min, length.out = len)
  if (log) 
    exp(lams)
  else lams
}

#' stability selection fit function for sparse inverse covariance using QUIC
#' @param x data matrix
#' @param y data matrix
#' @param q number of variables
#' @details This is a wrapper for QUIC to be used in stability selection. Pass it as
#' the fitfunction to stabsel
#' @seealso stabsel
#' @return list with selected variables and path.
#' @examples 
#' if (require(huge)) {
#'  set.seed(10010)
#'  dat.hubs <- huge.generator(n=N, d=40, graph="hub")
#'  pcutoff <- 0.75
#'  PFER <- 10
#'  s.hubs <- stabsel(x=dat.hubs$data, fitfun=stabs.quic, 
#'                   cutoff=pcutoff, PFER=PFER)
#'
#' }
stabs.quic <- function(x, y, q, ...)
{
  ## sort out a lambda path
  if (!requireNamespace("QUIC")) {
    stop("Package ", sQuote("QUIC"), " is required but not available")
  }
  empirical.cov <- cov(x)
  max.cov <- max(abs(empirical.cov[upper.tri(empirical.cov)]))
  lams <- getLamPath(max.cov, max.cov*0.05, len=40)
  est <- QUIC::QUIC(empirical.cov, rho=1, path=lams,msg=0)
  ut <- upper.tri(empirical.cov)
  qvals <- sapply(1:length(lams), function(idx){
    m <- est$X[,,idx]
    sum(m[ut] != 0)
  })
  
  ## Not sure if it is better to have more or less than q
  lamidx <- which.max(qvals >= q)
  ## Need to return the entire upper triangle - think about how to save
  ## ram later
  M <- est$X[,,lamidx][ut]
  selected <- (M != 0)
  s <- sapply(1:lamidx, function(idx){
    m <- est$X[,,idx][ut] != 0
    return(m)
  })
  colnames(s) <- as.character(1:ncol(s))
  return(list(selected=selected, path=s))
}
