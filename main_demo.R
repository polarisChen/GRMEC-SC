source("E:\\R_code\\01 GRMEC-SC\\GRMEC-SC\\constructSv.R")
source("E:\\R_code\\01 GRMEC-SC\\GRMEC-SC\\clustering_update.R")

library(Matrix)
library(parallel)
library(aricode)
library(Seurat)


#------------------------main function-----------------------#
lambda1s = 10^c(-4:4)
lambda2s = 10^c(-4:4)

vec = vector()
all_result=matrix(nrow = length(lambda1s)*length(lambda2s),ncol = 8)
colnames(all_result)=c('NMI_Last','ARI_Last','NMI_KLast','ARI_KLast','NMI_Max','ARI_Max','NMI_Kmax','ARI_Kmax')
for (i in 1:length(lambda1s)) {
  for (j in 1:length(lambda2s)) {
    lambda1=lambda1s[i]
    lambda2=lambda2s[j]
    vec[((i-1)*9+j)]=paste0( '_Lam1_', lambda1, '_Lam2_', lambda2)
  }
}
row.names(all_result)=vec
rm(lambda1,lambda2)

# npc_values = c(25, 50, 75, 100, 125)

for(i in 1:length(lambda1s)){
  for(j in 1:length(lambda2s)){
    lambda1 = lambda1s[i]
    lambda2 = lambda2s[j]
    
    try({
      split_data = list()
      rna = read.csv("E:\\R_code\\01 GRMEC-SC\\GRMEC-SC\\data\\RNA.CSV",row.name=1)
      rna = LogNormalize(rna)
      adt = read.csv("E:\\R_code\\01 GRMEC-SC\\GRMEC-SC\\data\\ATAC.CSV",row.name=1)
      adt = LogNormalize(adt)
      
      split_data[[1]] = t(rna)
      split_data[[2]] = t(adt)
      gt = read.csv("E:\\R_code\\01 GRMEC-SC\\GRMEC-SC\\data\\truth.csv", row.names=1, stringsAsFactors=F)
      gt = as.vector(gt$trueType)  
  
      views = 2
      
      ##* cluster number estimation
      K = length(unique(gt))
      message("## Estimated cluster number:", K, "\n")
      
      
      #* construct W & neighbors 
      clusts = read.csv("E:\\R_code\\01 GRMEC-SC\\GRMEC-SC\\data\\cllMix_sv.csv")
      clusts = clusts[, 2:dim(clusts)[2]]
      res = apply(clusts, MARGIN = 2, function(i){
        N = length(unique(i))
        return(N)
      })
      res
    
      S  = constructSv("E:\\R_code\\01 GRMEC-SC\\GRMEC-SC\\data\\cllMix_sv.csv")
      
      ## dimension reduction(npc) calculation
      npc = 50
      message("## Caluculated dim reduction npc:", npc, "\n")
      
      # main function
      iter=100
      res <- clustering_update(split_data, gt, K, npc, S, lambda1=lambda1, lambda2=lambda2, 
                               iteration=iter)
      message("## output the result of nmi and ari……")
      print(res)
      
      all_result[((i-1)*9+j),1]<-res[["nmi"]][iter,1]
      all_result[((i-1)*9+j),2]<-res[["ari"]][iter,1]
      all_result[((i-1)*9+j),3]<-res[["nmi"]][iter,2]
      all_result[((i-1)*9+j),4]<-res[["ari"]][iter,2]
      
      all_result[((i-1)*9+j),5]<-max(res[["nmi"]][,1])
      all_result[((i-1)*9+j),6]<-max(res[["ari"]][,1])
      all_result[((i-1)*9+j),7]<-max(res[["nmi"]][,2])
      all_result[((i-1)*9+j),8]<-max(res[["ari"]][,2])
      
      # save res
      # saveRDS(res, paste0("E:\\R_code\\01 GRMEC-SC\\GRMEC-SC\\result\\", 'cellMix', '_Lam1_', lambda1, '_Lam2_', lambda2, '.rds'))
    })
  }}  
saveRDS(all_result,'E:\\R_code\\01 GRMEC-SC\\GRMEC-SC\\result\\cellMix_res.rds')

