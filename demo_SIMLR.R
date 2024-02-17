# required external packages for SIMLR
library("Matrix")
library("parallel")
library(aricode)
# load the palettes for the plots
library(grDevices)
library(SIMLR)


# test SIMLR.R on example
cat("Performing analysis for SIMLR","\n")

demo_SIMLR <- function(data){
  
  s = Sys.time()
  
  gt = data$gt   #需要转置，输入m*n
  expr = t(as.matrix(data$expr))
  res_example = SIMLR(X=expr ,c=length(unique(gt)))
  
  clust = res_example$y$cluster
  nmi = NMI(clust, gt)
  ari = ARI(clust, gt)
  
  time = Sys.time() - s
  
  return(list(clust=clust, nmi=nmi, ari=ari, time=time))
  
}

