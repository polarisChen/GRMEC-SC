constructSv <- function(path_impSv){
  clusts = read.csv(path_impSv)
  clusts = clusts[, 2:dim(clusts)[2]]
  Sv = lapply(clusts, function(clust){
    N = length(clust)
    map = unique(clust)
    K = length(map)
    
    S = matrix(0, N, N)
    for(r in 1:K){
      id = which(clust == map[r])
      S[id, id] = 1
    }
    return(S)
  })
  return(Sv)
}