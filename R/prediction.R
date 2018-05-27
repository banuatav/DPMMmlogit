prediction = function(testSample,  thetaStar, K_t, N_k, n, J, M_d, maxMd, Dcon, Dcat, Dstar, hyperPar, m){
  y_p       = testSample$y
  X_p       = testSample$X
  Dummy_p   = testSample$Dummy
  n_p       = testSample$n_pred
  prediction_prob = matrix(0, nrow = n_p, ncol = 2 * J + 1)
  
  X_con_p = as.matrix(X_p[,1:Dcon])
  X_cat_p = as.matrix(X_p[, -(1:Dcon)])
  R_p     = cbind(rep(1, n_p), X_con_p, Dummy_p)
  
  ll_x     = matrix(0, n_p, (K_t+1))
  ll_xy    = array(0, c(n_p, J, K_t+1))
  ll_y_x   = array(0, c(n_p, J, K_t+1))
  
  for(k in 1:K_t){
    
    # Parameters for this cluster
    n_k           = N_k[k]
    mu_j_k        = as.matrix(thetaStar$mu[,k])
    Sigma_j_k     = matrix(thetaStar$Sigma[,k], nrow = Dcon)
    pi_j_k        = matrix(thetaStar$pi[,k], ncol = Dcat)
    beta_j_k      = as.matrix(thetaStar$beta[,k])
    beta_matrix   = matrix(beta_j_k, ncol = (J-1), nrow = length(beta_j_k)/(J-1), byrow = TRUE)
    
    # prediction prob
    # Xbeta        = cbind(R_p %*% beta_matrix, rep(0, n_p))
    # maxXbeta     = max(Xbeta)
    # relXbeta     = Xbeta - maxXbeta  #MASimpl matrix - vector
    # denom        = log(rowSums(exp(relXbeta)))
    # ll_y_x[,,k]  = relXbeta - denom  #MASimpl matrix - vector
    ll_y_x[,,k]  = prob_mlogit(beta_matrix, R_p)
    ll_xcon      = log_dmvnorm_arma_mult(X_con_p, matrix(mu_j_k, nrow = n_p , ncol = Dcon, byrow = TRUE), matrix(c(Sigma_j_k), nrow = n_p , ncol = Dcon^2, byrow = TRUE))
    
    ll_xcat      = rep(0, n_p)
    for (d in 1:Dcat){
      ll_xcat = ll_xcat + as.matrix(log(pi_j_k[cbind( X_cat_p[,d], rep(d, length(X_cat_p[,d])))]))
    }
    ll_x[, k]    = ll_xcon + ll_xcat + log(n_k) - log(thetaStar$alpha + n)  #MASimpl vector - const
    ll_xy[,,k]   = ll_y_x[,,k] + ll_x[ ,k]
  }
  
  ll_xm = matrix(0, n_p, m)
  ll_y_xm = array(0, c(n_p, J, m))
  
  for(k in 1:m){
    Sigma0 = riwish_rcpp(hyperPar$nu_Sigma, hyperPar$nu_Sigma*hyperPar$upsilon_Sigma*diag(1,Dcon)) #~IW(nu_Sigma, nu_Sigma * upsilon_Sigma*I)
    mu0    = rmvnormrcpp(1, hyperPar$mu_mu, hyperPar$Sigma_mu) #~N(mu_mu, sigma_mu)
    pi0_vec    = rep(NA, Dcat* maxMd)
    for (d in 1:Dcat){ pi0_vec[((d-1)*maxMd+1):((d-1)*maxMd+M_d[d])] = rdirichlet(n = 1, alpha = (rep(hyperPar$a0[d], M_d[d])/M_d[d])) } # ~Dir(a0/M_d) 
    pi0 = matrix(pi0_vec, ncol = Dcat)
    beta0 = rmvnormrcpp(1, hyperPar$betabar, hyperPar$B) # ~N(betabar, B)
    beta_matrix = matrix(beta0, ncol = (J-1), nrow = length(beta0)/(J-1), byrow = TRUE)
    
    ll_y_xm[,,k]  =  prob_mlogit(beta_matrix, R_p)
    ll_xcon        = log_dmvnorm_arma_mult(X_con_p, matrix(mu0, nrow = n_p , ncol = Dcon, byrow = TRUE), matrix(c(Sigma0), nrow = n_p , ncol = Dcon^2, byrow = TRUE))
    ll_xcat        = matrix(0, n_p, 1)
    for (d in 1:Dcat){
      ll_xcat = ll_xcat +as.matrix(log(pi0[cbind( X_cat_p[,d], rep(d, length(X_cat_p[,d])))]))
    }
    ll_xm[, k]    = ll_xcon + ll_xcat +  log(thetaStar$alpha) - log(thetaStar$alpha+ n) #MASimpl vector - const
    
    }
  
  ll_x[, K_t+1]    = log(rowSums(exp(ll_xm)/m))
  ll_y_x[,,K_t+1]  = log(t(colSums(exp(aperm(ll_y_xm, c(3,2,1))))/m))
  ll_xy[,,K_t+1]   = ll_y_x[,,K_t+1] + ll_x[ ,K_t+1]
  
  
  m_x                      = max(ll_x)
  # ll_xRel                  = ll_x - m_x
  # ll_xr                    = exp(ll_xRel)
  # ll_xr                    = ll_xr/rowSums(ll_xr)
  
  sum_xy                    = t(colSums(aperm(exp(ll_xy), c(3,2,1))))
  ll_y_x[,,1:K_t]             = ll_y_x[,,1:K_t] + log(N_k[1:K_t]) - log(thetaStar$alpha+ n)
  ll_y_x[,,K_t+1]             = ll_y_x[,,K_t+1]+ log(thetaStar$alpha) - log(thetaStar$alpha+ n)
  sum_y_x                     = t(colSums(aperm(exp(ll_y_x), c(3,2,1))))
  
  #sum_x                    = rowSums(exp(ll_x))
  sum_x                     = exp(m_x+log(rowSums(exp(ll_x - m_x))))
  #sum_xyNeal                = sum_x * sum_xyNeal
  
  prediction_prob[ , 1:J] = prediction_prob[ , 1:J] + sum_xy
  prediction_prob[ , J+1] = prediction_prob[ , J+1] + sum_x
  prediction_prob[ , (J+2):(2*J + 1)] = prediction_prob[ , (J+2):(2*J + 1)] + sum_y_x
  #predictionprob[ , (2*J+2):(3*J + 1)] = predictionprob[ , (2*J+2):(3*J + 1)] + sum_xyNeal
  
  ind = which(sum_x == 0)
  if(sum(sum_xy[-ind, ]/sum_x[-ind])!=n_p){
    warnings("Prediction probabilities do not sum to one.")
    if(sum(is.infinite(prediction_prob))!=0){warnings("Prediction probabilities contain infinite values")}
    if(sum(is.na(prediction_prob))!=0){warnings("Prediction probabilities contain NA")}
  }
  
  return(prediction_prob)
}
