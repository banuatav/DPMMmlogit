sampleCluster = function(thetaStar, clusterInfo, DAInfo, hyperPar){
  
  z_t   = clusterInfo$z   #J
  K_t   = clusterInfo$K   #
  N_k_t = clusterInfo$N_k #nj
  
  Dcon  = DAInfo$Dcon 
  Dcat  = DAInfo$Dcat
  Dstar = DAInfo$Dstar
  m     = DAInfo$m
  J     = DAInfo$J 
  n     = DAInfo$n
  M_d   = DAInfo$M_d 
  maxMd = DAInfo$maxMd 
  cumindex_xcat = seq(0,maxMd*Dcat-1, maxMd)
  UPDATE = TRUE
  
  # Initialize arrays for aux vars
  muAux    = matrix(NA, nrow = Dcon, ncol = m)
  SigmaAux = matrix(NA, nrow = Dcon*Dcon, ncol = m)
  piAux    = matrix(NA, nrow = Dcat*maxMd, ncol = m)
  betaAux  = matrix(NA, nrow = Dstar, ncol =  m)
  
  nu_Sigma      = hyperPar$nu_Sigma
  upsilon_Sigma = hyperPar$upsilon_Sigma
  #SigUps        = nu_Sigma*upsilon_Sigma*diag(1,Dcon)
  mu_mu         = hyperPar$mu_mu
  Sigma_mu      = hyperPar$Sigma_mu
  betabar       = hyperPar$betabar
  B             = hyperPar$B
  a0            = hyperPar$a0
  
  chol_B        = chol_rcpp(B)
  chol_Sigma_mu = chol_rcpp(Sigma_mu)
  chol_inv_SigUps =  1/sqrt(upsilon_Sigma*nu_Sigma) * diag(1,Dcon)
  
  for(i in 1:n){
    
    cluster_zi = z_t[i]
    K.min = length(unique(z_t[-i]))
    
    if (K_t == K.min){ #no cluster removed 
      # remove cluster assignment ind i
      z_t[i] = NA
      # Remove one from cluster of i 
      N_k_t[cluster_zi] = N_k_t[cluster_zi] - 1
      
      # muAux   = t(rmvnormrcpp( m, mu_mu, Sigma_mu))
      # betaAux = t(rmvnormrcpp(m, betabar, B))
      # for (d in 1:Dcat){piAux[((d-1)*maxMd+1):((d-1)*maxMd+M_d[d]), ] = t(rdirichlet(n = m, alpha = rep(a0[d], M_d[d])/M_d[d])) }      # ~Dir(a0) 
      # # Draw new values clusters k.min+1 to k.min+m
      # for (k in 1:m){SigmaAux[, k] = c(riwish_rcpp(nu_Sigma, SigUps))} #~IW(nu_Sigma, nu_Sigma * upsilon_Sigma*I)
      # 
      muAux   = t(rmvnormrcpp_chol( m, mu_mu, chol_Sigma_mu))
      betaAux = t(rmvnormrcpp_chol(m, betabar, chol_B))
      for (d in 1:Dcat){piAux[((d-1)*maxMd+1):((d-1)*maxMd+M_d[d]), ] = t(rdirichlet(n = m, alpha = rep(a0[d], M_d[d])/M_d[d])) }      # ~Dir(a0) 
      # Draw new values clusters k.min+1 to k.min+m
      for (k in 1:m){SigmaAux[, k] = c(riwish_rcpp_chol(nu_Sigma, chol_inv_SigUps))} #~IW(nu_Sigma, nu_Sigma * upsilon_Sigma*I)
      
    } else { 
      #cluster cluster_zi is removed: relabel clusters  and draw new values k.min+1 to k.min+m
      # relabel cluster assignment z_t. All cluster nr > cluster_zi 
      z_t[z_t>cluster_zi]        = z_t[z_t>cluster_zi] - 1
      # Cluster ind i now not available
      z_t[i]                      = NA
      # relabel number of obs in clusters 
      N_k_t[cluster_zi:(K.min+1)] = N_k_t[(cluster_zi+1):(K.min+2)] 
      
      # First aux cluster parameter is removed cluster's parameter
      SigmaAux[, 1] = thetaStar$SigmaStar[, cluster_zi]
      muAux[, 1]    = thetaStar$muStar[,cluster_zi]
      piAux[ , 1]   = thetaStar$piStar[, cluster_zi]
      betaAux[, 1]  = thetaStar$betaStar[,cluster_zi]
      
      # Adjust for changes in thetaStar: move all clusterparameters to match relabeled clusterassignments
      thetaStar$muStar[, cluster_zi:(K.min+1)]    = thetaStar$muStar[ , (cluster_zi+1):(K.min+2)] 
      thetaStar$SigmaStar[, cluster_zi:(K.min+1)] = thetaStar$SigmaStar[ , (cluster_zi+1):(K.min+2)] 
      thetaStar$piStar[, cluster_zi:(K.min+1)]    = thetaStar$piStar[ , (cluster_zi+1):(K.min+2)] 
      thetaStar$betaStar[, cluster_zi:(K.min+1)]  = thetaStar$betaStar[ , (cluster_zi+1):(K.min+2)] 
      
      # Draw new values clusters k.min+2 to k.min+m
      if(m>1){
        # muAux[,-1]    = t(rmvnormrcpp( m-1, mu_mu, Sigma_mu))
        # betaAux[,-1]  = t(rmvnormrcpp(m-1, betabar, B)) # ~N(betabar, B)
        # for (d in 1:Dcat){ piAux[((d-1)*maxMd+1):((d-1)*maxMd+M_d[d]), -1] = t(rdirichlet(n = m-1, alpha = rep(a0[d], M_d[d])/M_d[d])) } # ~Dir(a0) 
        # for (k in 2:m){SigmaAux[, k] = c(riwish_rcpp(nu_Sigma, SigUps))}
        # 
        muAux[,-1]    = t(rmvnormrcpp_chol(m-1, mu_mu, chol_Sigma_mu))
        betaAux[,-1]  = t(rmvnormrcpp_chol(m-1, betabar, chol_B))
        for (d in 1:Dcat){piAux[((d-1)*maxMd+1):((d-1)*maxMd+M_d[d]), -1] = t(rdirichlet(n = m-1, alpha = rep(a0[d], M_d[d])/M_d[d])) }      # ~Dir(a0) 
        # Draw new values clusters k.min+1 to k.min+m
        for (k in 2:m){SigmaAux[, k] = c(riwish_rcpp_chol(nu_Sigma, chol_inv_SigUps))} #~IW(nu_Sigma, nu_Sigma * upsilon_Sigma*I)
        
      }
      UPDATE = TRUE
    }
    
    # Sample probabilities for q1=existing clusters, q2=aux clusters
    q1 = rep(0, K.min)
    q2 = rep(0, m)
    
    # Compute loglikelihood given cluster for obs i
    y_i     = DAInfo$y[i]
    x_icon  = DAInfo$X[i, 1:Dcon]
    x_icat  = DAInfo$X[i, -(1:Dcon)]
    dum_i   = DAInfo$Dummy[i, ]
    
    if(UPDATE == TRUE){
      muS    = as.matrix(thetaStar$muStar[ ,1:K.min])
      SigmaS = as.matrix(thetaStar$SigmaStar[ ,1:K.min])
      piS    = as.matrix(thetaStar$piStar[ ,1:K.min])
      betaS  = as.matrix(thetaStar$betaStar[ ,1:K.min])
      UPDATE = FALSE
    }
    
    beta   = cbind(betaS, betaAux)
    r_i   = rbind(matrix(c(1, x_icon, dum_i), nrow = 1) %x% diag(1, J-1), rep(0,Dstar))
    Xbeta = r_i %*% beta
    maxXb = apply(Xbeta, 2, max)
    Xbeta = Xbeta - matrix(maxXb, nrow = J, ncol = K.min + m, byrow = TRUE) #MASimpl vector + num
    denom = log(colSums(exp(Xbeta)))
    
    ll_pij = Xbeta[y_i,] - denom 
    q1 = ll_pij[1:K.min]
    q2 = ll_pij[-(1:K.min)]
    
    index_xicat = cumindex_xcat + x_icat 
    if(K.min ==1){
      q1 = q1 + sum(log(piS[index_xicat,]))
    }else{
      q1 = q1 + colSums(log(piS[index_xicat,]))
    }
    if(m ==1){
      q2 = q2 + sum(log(piAux[index_xicat,]))
    }else{
      q2 = q2 + colSums(log(piAux[index_xicat,]))
    }
    
    q1 = q1 + c(log_dmvnorm_arma_mult(matrix(x_icon, nrow = K.min, ncol = Dcon, byrow = TRUE), t(muS), t(SigmaS)))
    q2 = q2 + c(log_dmvnorm_arma_mult(matrix(x_icon, nrow = m, ncol = Dcon, byrow = TRUE), t(muAux), t(SigmaAux)))
    
    # Adjust for probability 'in front' of likelihood
    q1 = q1 + log(N_k_t[1:K.min]) - log(n - 1 + thetaStar$alpha) #MASimpl vector + const
    q2 = q2 + log(thetaStar$alpha) - log(m) - log(n - 1 + thetaStar$alpha) #MASimpl vector + const
    
    # combine probabilities and 'un-log'
    q = c(q1, q2)
    maxq = max(q)
    qrel = q - maxq #added so that the numbers to be summed arent too big 
    q = exp(qrel) #exp(lli)/sum(exp(llj))=exp(lli^r +max)/sum(exp(llj^r +max)) = exp(lli^r)exp(max)/exp(max)sum(exp(llj^r)) =exp(lli^r)/sum(exp(llj^r)) 
    
    # calculate cumulative probabilities and sample from q
    q = q/sum(q)
    cumSum = cumsum(q)
    u = runif(n = 1, min = 0, max = 1)
    picked =  which(cumSum>u)
    picked = picked[1]
    
    # assign value z_i and its cluster parameters
    if (picked <= K.min){ # z_i belongs to an existing cluster
      z_t[i]        = picked 
      N_k_t[picked] = N_k_t[picked] + 1 
      K_t           = K.min
    } else{ # new cluster has opened
      K_t                         = K.min + 1
      z_t[i]                      = K_t 
      N_k_t[K_t]                  = 1
      
      picked_aux                  = picked - K.min
      
      # Adjust L for changes in thetaStar
      thetaStar$muStar[,K_t]    = muAux[,picked_aux]
      thetaStar$SigmaStar[,K_t] = SigmaAux[,picked_aux]
      thetaStar$piStar[,K_t]    = piAux[,picked_aux]
      thetaStar$betaStar[,K_t]  = betaAux[,picked_aux]
      UPDATE = TRUE
    }
  }
  
  # Adjust 
  clusterInfo$z   = z_t
  clusterInfo$N_k = N_k_t
  clusterInfo$K   = K_t
  
  return(list(thetaStar = thetaStar, clusterInfo = clusterInfo))
}