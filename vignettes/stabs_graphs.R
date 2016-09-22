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
## this is the wrapper function for QUIC
quic.graphical_model

## ---- StabsRun ----
pcutoff <- 0.75
PFER <- 10
s.hubs <- stabsel(x=dat.hubs$data, fitfun=quic.graphical_model, 
                  cutoff=pcutoff, PFER=PFER)
s.cluster <- stabsel(x=dat.cluster$data, fitfun=quic.graphical_model, 
                  cutoff=pcutoff, PFER=PFER)
s.rand <- stabsel(x=dat.rand$data, fitfun=quic.graphical_model, 
                  cutoff=pcutoff, PFER=PFER)

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
