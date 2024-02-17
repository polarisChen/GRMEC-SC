library(cidr)
library(aricode)


demo_CIDR <- function(data){   #??Ҫת?ã?????m*n
  
  s = Sys.time()
  
  par(ask=FALSE)
  
  gt = data$gt
  colnames(data$expr)=c(1:dim(data$expr)[2])
  if(class(data$expr) == "dgCMatrix"){
    data$expr = as.matrix(data$expr)
  }
  data$expr<-na.omit(data$expr)
  #*data = scDataConstructor(data)
  data1 = scDataConstructor(t(data$expr))
  data1 = determineDropoutCandidates(data1)
  data1 = wThreshold(data1)
  data1@wThreshold <- 7    # 修改权重阈值
  data1 = scDissim(data1)
  #*
  data1 = scPCA(data1,plotPC = F)
  data1 = nPC(data1)
  nCluster(data1)
  data1 = scCluster(data1, nCluster = 4)  # 细胞类型
  clust = data1@clusters
  
  nmi = NMI(clust, gt)
  ari = ARI(clust, gt)
  
  
  time = Sys.time() - s
  
  result = list(clust=clust, nmi=nmi, ari=ari, time=time)
  
  return(list(clust=clust, nmi=nmi, ari=ari, time=time))
}
