sampleTheta = function(thetaStar, clusterInfo , DAInfo, hyperPar){
  
  Dcon      = DAInfo$Dcon 
  Dcat      = DAInfo$Dcat 
  Dstar     = DAInfo$Dstar
  J         = DAInfo$J
  M_d       = DAInfo$M_d
  maxMd       = max(M_d)
  
  z = clusterInfo$z 
  K = clusterInfo$K
  
  mu_mu    = hyperPar$mu_mu
  upsilon_mu = hyperPar$upsilon_mu
  nu_mu = hyperPar$nu_mu
  Sigma_mu = hyperPar$Sigma_mu
  
  upsilon_Sigma = hyperPar$upsilon_Sigma
  nu_Sigma =  hyperPar$nu_Sigma
  
  a0       = hyperPar$a0 
  betabar  = hyperPar$betabar 
  B        = hyperPar$B 
  s        = DAInfo$s
  
  for( k in 1:K){
    # Count number of times cluster has been open
    DAInfo$ClusterOpen[k] = DAInfo$ClusterOpen[k] + 1
    
    ## Select Data
    y_k  = DAInfo$y[z==k]
    Nk   = length(y_k)
    
    if ( Nk == 1 ) {
      X_k         = t(as.matrix(DAInfo$X[z==k, ]))
      X_k_con     = t(as.matrix(X_k[,1:Dcon]))
      X_k_cat     = t(as.matrix(X_k[,-(1:Dcon)]))
      Dummy_k     = t(as.matrix(DAInfo$Dummy[z == k, ]))
      Dy_k        = t(as.matrix(DAInfo$Dy[z == k, ]))
    } else {
      X_k  = DAInfo$X[z==k, ] 
      X_k_con     = as.matrix(X_k[,1:Dcon])
      X_k_cat     = as.matrix(X_k[,-(1:Dcon)])
      Dummy_k     = DAInfo$Dummy[z == k, ]
      Dy_k        = DAInfo$Dy[z == k, ]
    }
    
    mu_k      =  as.matrix(thetaStar$muStar[, k])
    iota      =  as.matrix(rep(1,Nk))
    ###################### 2.1.1 Sample sigma ######################
    # update parameters Sigma_k ~ IW(df_draw, Tau_draw): df_draw = Nk + nu0 +1, Tau_draw =upsilon0*I_Dcon + c*SumSqmu + SumSqX
    diff_X      = X_k_con - tcrossprod(iota, mu_k)
    S_x         = crossprod(diff_X)
    df_draw     = nu_Sigma + Nk
    Tau_draw    = nu_Sigma * upsilon_Sigma*diag(1, Dcon) + S_x
    
    # new draw Sigma_k
    Sigma_new = riwish_rcpp(df_draw, Tau_draw)
    
    ###################### 2.1.2 Sample mu ######################
    # update parameters mu_k ~ N((cmu0 + sumX_over_i/(c+nk), (c+nk)^-1 Sigma_k)
    inv_mu_Sigma = chol2inv(chol(Sigma_mu))
    inv_Sigma_k  = chol2inv(chol(Sigma_new))
    Sigma_draw  = chol2inv(chol(inv_mu_Sigma + Nk *inv_Sigma_k))
    sumX        = as.matrix(colSums(X_k_con))
    mu_draw     = Sigma_draw %*% (inv_mu_Sigma %*% as.matrix(mu_mu) + inv_Sigma_k %*% sumX)
    
    # new draw mu_k
    mu_new = rmvnormrcpp(1, mu_draw, Sigma_draw)
    
    ###################### 2.1.3 Sample pi_d ######################
    pi_new = rep(NA, maxMd*Dcat)
    
    # update parameters and draw pi_d,k ~ Dir(a0_d/M_d + sum(I[x_i,d == 1], ... , a0_d/M_d + sum(I[x_i,d == Md]))
    for (d in 1:Dcat){
      a_draw = rep(0,M_d[d])
      a = a0[d]/ M_d[d]
      for (m  in 1:M_d[d]){a_draw[m] = a + sum(X_k_cat[ , d]==m)}
      pi_new[((d-1)*maxMd+1):((d-1)*maxMd+M_d[d])] = rdirichlet(n = 1, alpha = a_draw)
    }
    
    # Update thetaStar
    thetaStar$muStar[ , k]    = mu_new
    thetaStar$SigmaStar[ , k] = c(Sigma_new)
    thetaStar$piStar[ , k]    = pi_new
    
    ###################### 2.2 Sample Theta*^y=beta ######################
    betacurrent       = thetaStar$betaStar[ , k]
    R_k               = cbind(rep(1, Nk), X_k_con, Dummy_k)
    
    betaMat = matrix(0, nrow = Dstar/(J-1), ncol = J)
    betaMat[,-J] = matrix(betacurrent, ncol = J-1, byrow = TRUE)
    B_N = (s^2)*invHess(B, y_k, R_k, betaMat, J)
    
    beta_candidate = betacurrent + as.numeric(rmvnormrcpp(1, rep(0, Dstar), B_N))
    
    llik_cur = getLogLikMnl(betacurrent, Dy_k, R_k, J)
    llik_can = getLogLikMnl(beta_candidate, Dy_k, R_k, J)
    
    d_log_l     = llik_can - llik_cur
    d_log_prior = log_dmvnorm_arma_sing(beta_candidate, betabar, B) - log_dmvnorm_arma_sing(betacurrent, betabar, B)
    
    log_alpha =  d_log_l  +  d_log_prior 
    alpha_min = min(exp(log_alpha), 1)
    u = runif(n = 1, min = 0, max = 1)
    if(u<=alpha_min){
      thetaStar$betaStar[ , k] = beta_candidate
      DAInfo$AcceptRate[k] = DAInfo$AcceptRate[k] + 1
    } 
    
  }
  return(list(thetaStar=thetaStar, AcceptRate = DAInfo$AcceptRate, ClusterOpen = DAInfo$ClusterOpen ))
}