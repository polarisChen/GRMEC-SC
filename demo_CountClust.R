library(maptpx)
library(aricode)
library(CountClust)

demo_CountClust <- function(data, K){

  #* count transform for imp_X
  #data$expr = round(2^(data$expr) - 1)
  
  s = Sys.time()

  gene_names = colnames(data$expr)
  meta_data = rownames(data$expr)
  gt = data$gt
  
  
  # return matrices omega--cluster indicator & theta--basis matrix
  if(class(data$expr) == "dgCMatrix" | class(data$expr) == "data.frame"){
    data$expr = as.matrix(data$expr)
  }
  
  # remove zero gene
  id = which(apply(data$expr, 2, function(i){all(i==0)}))
  result = FitGoM(data$expr, K=K)
  clust = apply(result$fit$omega, 1, which.max)
  names(clust) = meta_data
  nmi = NMI(clust, gt)
  ari = ARI(clust, gt)
  
  if(FALSE){
    omega = result$fit$omega
    annotation = data.frame(
      sample_id = paste0("X", c(1:NROW(omega))),
      tissue_label = factor(rownames(omega),
                            levels = rev( c("zy", "early2cell",
                                            "mid2cell", "late2cell",
                                            "4cell", "8cell", "16cell",
                                            "earlyblast","midblast",
                                            "lateblast") ) ) )
    
    
    StructurePie(t(data$expr), input_type="apply_pca",
                 use_voom=F, omega=omega, xlab="PCA1",
                 ylab="PCA2",
                 main="structure k=3 pie on PCA")
    
    
    theta_mat = result$fit$theta
    top_features = ExtractTopFeatures(theta_mat, top_features=100,
                                      method="poisson", options="min")
    gene_list = do.call(rbind, lapply(1:dim(top_features)[1],
                                      function(x) gene_names[top_features[x, ]]))
  }

  time = Sys.time() - s
  return(list(clust = clust, nmi = nmi, ari = ari, time = time))
  
}

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
