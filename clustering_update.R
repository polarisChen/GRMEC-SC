#####################################################
### update function of PDGE_EC
### cost function -- \sum wv||Xv - W Vv|| + ??1tr(W' LH W) + ??2 \sum ??vtr(H' Lsv H)



clustering_update <- function(X, gt, K, npc, S, lambda1, lambda2, iteration=500, 
                              thre=1e-4, rep=1){
  s = Sys.time()
  cat("## clustering vars updating...\n")
  N = dim(X[[1]])[1]  #????????ͬһ??N
  P = dim(X[[1]])[2]
  Mx = length(X)
  Ms = length(S)

  # init vars
  set.seed(rep)   #*15\ produce reproductive result
  H = eigen(diag(colSums(S[[1]])) - S[[1]])$vector
  H = H[, (N-K):(N-1)]
  H[H<=0] = 1e-10
  W = matrix(runif(N * npc), N, npc) 
  V = lapply(1:Mx, function(i){    
    pi = dim(X[[i]])[2]
    Vi = matrix(runif(npc * pi), npc, pi)
    return(Vi)
  })
  Sen = matrix(runif(N * N), N, N)
  L = lapply(S, function(i){
    res = diag(rowSums(i)) - i
    return(res)
  })
  
  
  # initialization
  J_set = NULL
  J_DR = matrix(0, iteration, Mx)
  J_MV = matrix(0, iteration, Ms)
  J_LE = NULL
  J_HE = NULL
  wv = matrix(1/Mx, iteration, Mx)
  alpha = matrix(1/Ms, iteration, Ms)
  par1 = NULL
  par2 = NULL
  par_sigma = NULL
  par_beta = NULL
  nmi = matrix(0, iteration, 2)
  ari = matrix(0, iteration, 2)

  
  for(iter in 1:iteration){
    cat(paste0("### Updating ", iter, "th round of ", iteration, "\n"))
    cat(paste0(strrep("#", round(30*iter/iteration)), " ", 
               100*signif(iter/iteration, digits=2), "%\n"))

    #-------------------- update Wv and alphav ---------------------#
    cat('### updating wv and alphav...\n')
    J_DR[iter, ] = t(sapply(1:Mx, function(v){
      cost = norm(X[[v]] - W %*% V[[v]], 'F')^2
      return(cost)
    }))
    wv[iter, ] = t(sapply(1:Mx, function(v){
      res = 0.5 / sqrt(J_DR[iter, v])
      return(res)
    }))
    # 4\unitize weight
    wv[iter, ] = wv[iter, ] / sum(wv[iter, ])
    
    J_MV[iter, ] = t(sapply(1:Ms, function(v){
      J_MVv = norm(Sen - S[[v]], "F")^2
      return(J_MVv)
    }))
    alpha[iter, ] = t(sapply(1:Ms, function(v){
      res = 0.5 / sqrt(J_MV[iter, v])
      return(res)
    }))
    #alpha[iter, ] = alpha[iter, ] / sum(alpha[iter, ])
    
    #-------------------- update W ---------------------#
    cat("### updating W...\n")
    SH = H %*% t(H)
    LH = diag(colSums(SH)) - SH
    DH = diag(colSums(SH))
    mvpr = 0
    mvnw = 0
    for(v in 1:Mx){
      mvpr = mvpr + wv[iter, v] * V[[v]] %*% t(V[[v]])
      mvnw = mvnw + wv[iter, v] * X[[v]] %*% t(V[[v]])
    }
    fir_ord_g = mvnw + lambda1 * SH %*% W
    sec_ord_g = lambda1 * DH %*% W + W %*% mvpr
    W = W * (fir_ord_g / sec_ord_g)
    rm(mvpr, mvnw, SH)

    
    # 5\ *normalize
    norms = rowSums(W)
    norms[norms == 0] = 1e-10
    W = W / matrix(norms, N, dim(W)[2])
    
    
    #-------------------- update Vv ---------------------#
    cat("### updating V...\n")
    V = lapply(1:Mx, function(v){
      Vv = V[[v]] * ((t(W) %*% X[[v]]) / (t(W) %*% W %*% V[[v]]))
      return(Vv)
    })
    
    
    #-------------------- update Sen ---------------------#
    cat("### updating Sen...\n")
    dH = rowSums(H^2)
    dH = matrix(dH, nrow=N, ncol=N)
    dH = dH + t(dH) - 2*H %*% t(H)
    dH = dH - diag(diag(dH)) + diag(1e-10, N) #* remove effect from self sim
    
    temp = lapply(1:Ms, function(v){
      res = alpha[iter, v] * S[[v]]
      return(res)
    })
    
    #**
    Sen = Sen * (2 * Reduce(f="+", temp)) /
      (2 * sum(alpha[iter, ]) * Sen + lambda2 / 2 * dH + 1e-10)
    Sen[which(Sen == 0)] = 1e-10
    
    
    #-------------------- update H ---------------------#
    cat("### updating H...\n")
    d = rowSums(W^2)
    d = matrix(d, nrow=N, ncol=N)
    d = d + t(d) - 2*W %*% t(W)
    
    H = H * ((lambda2 * Sen %*% H) / ((lambda2 * diag(rowSums(Sen)) + lambda1 / 2 * d) %*% H))
    #* 5\row based H normalize
    norms = rowSums(H)    
    norms[norms == 0] = 1e-10
    H = H / matrix(norms, N, K)

    
    #-------------------- update paras ---------------------#
    #** 6\update parameters    
    if(FALSE){#iter %% 10 == 0){
      He = sqrt(eigen(H %*% t(H))$values)[1]
      Ve = c(sqrt((eigen(V[[1]] %*% t(V[[1]]))$values)[1]), sqrt((eigen(V[[2]] %*% t(V[[2]]))$values)[1]),
             sqrt((eigen(V[[3]] %*% t(V[[3]]))$values)[1]))
      lambda1 = sum(wv[iter, ] * Xe) / (He + sum(wv[iter, ] * Ve))
      
      de = eigen(d)$values[1]
      lambda2 = (lambda1 * de) / (2 * Se)
      par1[iter] = lambda1
      par2[iter] = lambda2
    }
    
    #-------------------- update cost ---------------------#
    J_HE[iter] = sum(diag(t(W) %*% LH %*% W))
    J_LE[iter] = sum(diag(t(H) %*% (diag(rowSums(Sen))-Sen) %*% H))
    J_set[iter] = sum(wv[iter, ] * J_DR[iter, ]) + lambda1 * J_HE[iter] + sum(alpha[iter, ] * J_MV[iter, ]) + lambda2 * J_LE[iter]
    cat("### Current cost:", J_set[iter], "\n")
    
    
    if(iter > 1){
      var_J = abs(J_set[iter] - J_set[iter-1])
      
      # convergence check
      if(FALSE){#var_J <= thre){
        cluster = sapply(1:N, function(i){
          which.max(H[i, ])
        })
        
        res = list(W = W, V = V, H = H, cluster = cluster, Dc = diag(Dc), Dg = diag(Dg),
                   J = J_set, J_DR = J_DR, J_HE = J_HE, J_LE = J_LE, nmi = nmi, ari = ari,
                   lambda1 = par1, lambda2 = par2, wv = wv, dw = d)   #* lambda save

        return(res)
      }
    }
    
    
    # recording nmi
    cluster = sapply(1:N, function(i){
      which.max(H[i, ])
    })
    nmiH_max = NMI(cluster, gt)
    ariH_max = ARI(cluster, gt)
    nmi[iter, 1] = nmiH_max
    ari[iter, 1] = ariH_max
    tryCatch(
      {
        cluster_k = kmeans(H, K, nstart=20, iter.max = 200)$cluster
        nmiH_kmean = NMI(cluster_k, gt)
        ariH_kmean = ARI(cluster_k, gt)
        cat(paste0("## NMI:   nmiH_max: ", nmiH_max, "    nmiH_kmean: ", nmiH_kmean, "\n"))
        cat(paste0("## ARI:   ariH_max: ", ariH_max, "    ariH_kmean: ", ariH_kmean, "\n"))
        nmi[iter, 1:2] = c(nmiH_max, nmiH_kmean)
        ari[iter, 1:2] = c(ariH_max, ariH_kmean)
      }, error=function(e){
        cat(paste0("## NMI:   nmiH_max: ", nmiH_max, "\n"))
        cat(paste0("## ARI:   ariH_max: ", ariH_max, "\n"))
      }
    )

  }
  
  # saving result
  cat("# MGGE iteration complete!\n")
  time = Sys.time() - s
  message("## Consume", time, "seconds.\n")

  
  res = list(gt = gt, W = W, V = V, H = H, cluster = cluster, S = S, Sen = Sen, 
             J = J_set, J_DR = J_DR, J_HE = J_HE, J_LE = J_LE, J_MV = J_MV, nmi = nmi, ari=ari,
             lambda1 = lambda1, lambda2 = lambda2, wv = wv, alpha = alpha, dw = d)   #* lambda save

  return(res)
}


soft <- function(x, theta){
  if(max(theta == 0)){
    res = x
    return(x)
  }else{
    res = abs(x) - theta
    res[which(as.matrix(res<0))] = 0
    res = sign(x) * res
    return(res)
  }
}
