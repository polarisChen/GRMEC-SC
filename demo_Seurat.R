library(Seurat)
library(aricode)
library(Matrix)


demo_Seurat <- function(data, resolution){
  
  s = Sys.time()
  gt = data$gt
  cell_name = rownames(data$expr)
  gene_name = colnames(data$expr)
  N = length(cell_name)
  G = length(gene_name) 
  data = Matrix(log2(t(data$expr) + 1), sparse = TRUE)
  
  data = CreateSeuratObject(counts = data, min.cells = 3, min.features = 1)
  
  # preprocessing
#  data[["percent.mt"]] = PercentageFeatureSet(data, pattern="^MT-")
  data = NormalizeData(object = data, normalization.method = "LogNormalize", 
                       scale.factor = 1e4, display.progress = FALSE)
  
  # find hv_gene
  # if(G < 300){
  #   data <- FindVariableFeatures(data, selection.method = "vst", nfeatures =G)
  # }else{
  #   data = FindVariableFeatures(data, selection.method = "vst", nfeatures =2000)#2000
  # }
  
  if(G > 300){
    data = FindVariableFeatures(data, selection.method = "vst", nfeatures =2000)#2000
  }else{
    data = FindVariableFeatures(data)
  }
  
  # scaling the data
  data <- ScaleData(object = data, display.progress = FALSE)
  
  # dimensional reduction
  # if(N < 300){
  #   data = RunPCA(data, features=VariableFeatures(object=data),
  #                 npcs=N-1)
  # }else{
  #   data = RunPCA(data, features=VariableFeatures(object=data))
  # }
  data <- RunPCA(object = data,  verbose = FALSE)
  # Estimate the num of cluster
  # data = JackStraw(data, num.replicate=100)
  # data = ScoreJackStraw(data, dims=1:20)

  # supervise the num of cluster with alternative methods
  # JackStrawPlot(data, dims=1:15)
  # ElbowPlot(data)
  

  # clustering
  if(G < 300){
    data = FindNeighbors(data, dims=1:6)
    res = lapply(resolution, function(i){
      data = FindClusters(data, resolution=i)
      clust = Idents(data)
      nmi = NMI(clust, gt)
      ari = ARI(clust, gt)
      return(list(clust=clust, nmi=nmi, ari=ari))
    })
  }else{
    data = FindNeighbors(data, dims=1:50)
    res = lapply(resolution, function(i){
      data = FindClusters(data, resolution=i)
      clust = Idents(data)
      nmi = NMI(clust, gt)
      ari = ARI(clust, gt)
      return(list(clust=clust, nmi=nmi, ari=ari))
    })
  }
  
  return(res)
}