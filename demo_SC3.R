library(SC3)
library(aricode)
library(SingleCellExperiment)


demo_SC3 <- function(data){
  
  s = Sys.time()
  
  gt = data$gt
  cell_name = rownames(data$expr)
  gene_name = colnames(data$expr)
  N = length(cell_name)
  G = length(gene_name)
  sce = SingleCellExperiment(
    assays = list(
      #counts = t(as.matrix(data$expr)),
      counts = t(as.matrix(data$expr)),
      logcounts = t(as.matrix(data$expr)+1)
      #logcounts = t(log2(as.matrix(data$expr) + 1))
    ),
    colData = data$gt
  )
  rowData(sce)$feature_symbol <- rownames(sce)
  
  k = length(unique(gt))
  h=k+1
  # sce <- sc3_prepare(sce, gene_filter = FALSE)
  # #sce = sc3_estimate_k(sce)
  # metadata(sce)$sc3$k_estimation = k
  # #metadata(sce)$sc3$k_estimation = min(metadata(sce)$sc3$k_estimation, 15)
  # #k = metadata(sce)$sc3$k_estimation
  # sce = sc3_calc_dists(sce)
  # sce = sc3_calc_transfs(sce)
  # sce = sc3_kmeans(sce, k)
  # sce = sc3_calc_consens(sce)
  
  if(N < 1999){
    if(G > 30){
      sceset <- sc3(sce, ks = k:h, biology = TRUE,gene_filter = T,pct_dropout_max = 98)#2000
      # sce <- sc3(sce,  ks = k:k+1, biology = TRUE, svm_num_cells = 50, n_cores = 1)#2000
      # sce <- sc3_run_svm(sce, ks = k)
    }else{
      sceset <- sc3(sce, ks = k:h, biology = TRUE,gene_filter=F)
    }
  }else{
    if(G > 30){
      sceset <- sc3(sce,  ks = k:h, biology = TRUE, svm_num_cells = 50)#2000
      sceset <- sc3_run_svm(sceset, ks = k)
    }else{
      sceset <- sc3(sce, ks = k:h, biology = TRUE, svm_num_cells = 50,gene_filter=F)
      sceset <- sc3_run_svm(sceset, ks = k)
    }
  }
  
  clust = colData(sceset)[, 2]
  
  nmi = NMI(clust, gt)
  ari = ARI(clust, gt)
  
  time = Sys.time() - s
  
  return(list(clust=clust, nmi=nmi, ari=ari, time=time))
  
}