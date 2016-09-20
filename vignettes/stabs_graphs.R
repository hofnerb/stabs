## ---- Setup ----
library(stabs)
if (!require(huge)) {
  stop("Need package huge for generating test data")
}
if (!require(QUIC)) {
  stop("Need package QUIC")
}
if (!require(igraph)) {
  stop("Need package igraph")
}
N <- 200
set.seed(10010)
dat.hubs <- huge.generator(n=N, d=40, graph="hub")
set.seed(10010)
dat.cluster <- huge.generator(n=N, d=40, graph="cluster")
set.seed(10010)
dat.rand <- huge.generator(n=N, d=40, graph="random")

## ---- PlotHubs ----
plot(dat.hubs)

## ---- PlotClust ----
plot(dat.cluster)

## ---- PlotRand ----
plot(dat.rand)

## ---- StabsFunction ----
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

## ---- StabsRun ----
pcutoff <- 0.75
PFER <- 10
s.hubs <- stabsel(x=dat.hubs$data, y=dat.hubs$data, fitfun=stabs.quic, 
                  cutoff=pcutoff, PFER=PFER, B=200, graphical=TRUE)
s.cluster <- stabsel(x=dat.cluster$data, y=dat.cluster$data, fitfun=stabs.quic, 
                  cutoff=pcutoff, PFER=PFER, B=200, graphical=TRUE)
s.rand <- stabsel(x=dat.rand$data, y=dat.rand$data, fitfun=stabs.quic, 
                  cutoff=pcutoff, PFER=PFER, B=200, graphical=TRUE)

## ---- StabsPlot ----
p1 <- function(stabsout, orig)
{
  ## display comparison of original graph and stabs estimation
  j<- orig$omega * 0
  orig.graph <- graph.adjacency(orig$theta != 0, mode="max", diag=FALSE)
  ut <- upper.tri(j)
  j[ut][stabsout$selected] <- 1
  stabs.graph <- graph.adjacency(j!=0, mode="max", diag=FALSE)
  layout <- layout.fruchterman.reingold(orig.graph)
  par(mfrow=c(1,2))
  plot(orig.graph, layout=layout, edge.color="gray50", vertex.color='red', 
       main="Real graph",  vertex.size = 3, vertex.label = NA)
  plot(stabs.graph, layout=layout, edge.color="gray50", vertex.color='red', 
       main="Stabs estimated graph", vertex.size = 3, vertex.label = NA)
}
## ---- StabsPlotHubs ----
p1(s.hubs, dat.hubs)
## ---- StabsPlotCluster ----
p1(s.cluster, dat.cluster)
## ---- StabsPlotRand ----
p1(s.rand, dat.rand)
