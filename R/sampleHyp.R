sampleHyp =function(parHP, thetaStar, hyperPar, K, M_d, maxMd, Dcon, Dcat, Dstar){
  ## Samples the hyperparamters of continuous covariates which are all given uniform priors. 
  # The conditional posteriors are truncated at the lower and upperbound of the hyperpar distribution
  # Obtain information relevant for sampling
  
  mu_mu      = hyperPar$mu_mu
  Sigma_mu   = hyperPar$Sigma_mu
  upsilon_mu = hyperPar$upsilon_mu
  nu_mu      = hyperPar$nu_mu
  mu_mu_par  = parHP$mu_mu_par
  upsilon_mu_par = parHP$upsilon_mu_par
  
  nu_Sigma      = hyperPar$nu_Sigma
  upsilon_Sigma = hyperPar$upsilon_Sigma
  upsilon_Sigma_par = parHP$upsilon_Sigma_par
  
  a0_par = parHP$a0_par
  
  betabar     = hyperPar$betabar
  upsilonbeta = hyperPar$upsilonbeta
  nubeta      = hyperPar$nubeta
  betabar_par = parHP$betabar_par
  upsilonbeta_par = parHP$upsilonbeta_par
  
  muS    = as.matrix(thetaStar$muStar[,1:K])
  SigmaS = as.matrix(thetaStar$SigmaStar[,1:K])
  betaS  = as.matrix(thetaStar$betaStar[,1:K])
  piS    = as.matrix(thetaStar$piStar[,1:K])
  
  # Sample Sigma_mu from IW(nu_mu + K, nu_mu*upsilon_mu*I + S_mu), S_mu = sum (mu - mu_mu)(mu - mu_mu)^T
  df_SigMu_update  = nu_mu + K
  mu_dif = muS - matrix(mu_mu, nrow = Dcon, ncol = K)
  S_mu = tcrossprod(mu_dif)
  Tau_SigMu_update = nu_mu * upsilon_mu * diag(1, Dcon) + S_mu
  Sigma_mu_new = riwish_rcpp(df_SigMu_update, Tau_SigMu_update)
  
  mean_mu_mu_update = rowSums(muS)/K
  Sigma_mu_mu_update = Sigma_mu_new/K
  
  # Sample mu_mu from inverse cdf method on truncated distr
  mu_mu_new = mu_mu
  index = 1:Dcon
  for (d in 1:Dcon){
    index_notd = index[-d]
    xnotd = mu_mu_new[index_notd]
    
    mu_mu_d = mean_mu_mu_update[d]
    mu_mu_notd = mean_mu_mu_update[index_notd]
    
    cov_dd = Sigma_mu_mu_update[d,d]
    cov_notdnotd = Sigma_mu_mu_update[index_notd, index_notd]
    cov_notdd = Sigma_mu_mu_update[cbind(rep(d, length(index_notd)), index_notd)]
    
    A = crossprod(as.matrix(cov_notdd), chol2inv(chol(cov_notdnotd)))
    con_meand = mu_mu_d + A %*% as.matrix(xnotd - mu_mu_notd) 
    con_sigmad = cov_dd - A %*% as.matrix(cov_notdd)
    
    u = runif(1, 0, 1)
    lb = mu_mu_par[d, 1]
    ub = mu_mu_par[d, 2]
    Fz = pnorm(lb, con_meand, con_sigmad) + u*(pnorm(ub, con_meand, con_sigmad) - pnorm(lb, con_meand, con_sigmad))
    
    if(Fz != 1 & Fz != 0)
    {
      mu_mu_new[d] = qnorm(Fz, con_meand, con_sigmad)
    } else {
      mu_mu_new[d] = parHP$mu_mu_par[d, Fz+1]
    }
  }
  
  #  Update parameters cond poster upsilon_mu
  shape_upsilon_mu = 1 + nu_mu*Dcon/2
  rate_upsilon_mu = nu_mu*(tr(chol2inv(chol(Sigma_mu_new))))/2
  
  #  Sample upsilon_mu from inverse cdf method on truncated distr
  u = runif(1, 0, 1)
  lb = parHP$upsilon_mu_par[1]
  ub = parHP$upsilon_mu_par[2]
  Fz = pgamma(lb, shape = shape_upsilon_mu, rate = rate_upsilon_mu) + u*(pgamma(ub, shape = shape_upsilon_mu, rate = rate_upsilon_mu) - pgamma(lb,  shape = shape_upsilon_mu, rate = rate_upsilon_mu))
  if(Fz != 1 & Fz != 0)
  {
    upsilon_mu_new = qgamma(Fz,  shape = shape_upsilon_mu, rate = rate_upsilon_mu)
    
  } else {
    upsilon_mu_new = parHP$upsilon_mu_par[Fz+1]
  }
  
  
  trace_sumInv_SigmaS = 0
  for (k in 1:K){
    trace_sumInv_SigmaS = trace_sumInv_SigmaS + tr(chol2inv(chol(matrix(SigmaS[,k], nrow =2))))
  }
  
  #  Update parameters cond poster upsilon_Sigma
  shape_upsilon_Sigma = 1 + nu_Sigma*Dcon*K/2
  rate_upsilon_Sigma = nu_Sigma*(trace_sumInv_SigmaS)/2
  
  #  Sample upsilonbeta from inverse cdf method on truncated distr
  u = runif(1, 0, 1)
  lb = upsilon_Sigma_par[1]
  ub = upsilon_Sigma_par[2]
  Fz = pgamma(lb, shape = shape_upsilon_Sigma, rate = rate_upsilon_Sigma) + u*(pgamma(ub, shape = shape_upsilon_Sigma, rate = rate_upsilon_Sigma) - pgamma(lb,  shape = shape_upsilon_Sigma, rate = rate_upsilon_Sigma))
  if(Fz != 1 & Fz != 0)
  {
    upsilon_Sigma_new = qgamma(Fz,  shape = shape_upsilon_Sigma, rate = rate_upsilon_Sigma)
    
  } else {
    upsilon_Sigma_new = parHP$upsilon_Sigma_par[Fz+1]
  }
  
  # Sample new a0 from inverse cdf method on truncated distr
  a0_new = rep(0, Dcat)
  Fz = rep(NA, Dcat)
  
  for (d in 1:Dcat){
    sum_logpi_d = sum(log(piS[((d-1)*maxMd +1):((d-1)*maxMd +M_d[d]),]))
    rate_a0_d = -sum_logpi_d/M_d[d]
    
    u = runif(1, 0, 1)
    lb = parHP$a0_par[d, 1]
    ub = parHP$a0_par[d, 2]
    Fz[d] = pexp(lb, rate = rate_a0_d) + u*(pexp(ub, rate = rate_a0_d) - pexp(lb, rate = rate_a0_d))
    a0_new[d] = qexp(Fz[d], rate = rate_a0_d)
  }
  
  index_up = Fz == 1
  index_low = Fz == 0
  a0_new[index_up] = a0_par[index_up, 2]
  a0_new[index_low] = a0_par[index_low, 1]
  
  diff_beta = betaS - matrix(betabar, nrow = Dstar, ncol = K)
  sumbeta   = rowSums(betaS)
  sumSbeta  = tcrossprod(diff_beta)
  
  # Update parameters cond poster B
  nubeta_update = nubeta + K
  Sbeta_update = nubeta * upsilonbeta * diag(1,Dstar) + sumSbeta
  B_new = riwish_rcpp(nubeta_update, Sbeta_update)
  
  #  Update parameters cond poster betabar
  mean_betabar_update = sumbeta/K
  Sigmabetabar_update = B_new/K
  
  # Sample betabar from inverse cdf method on truncated distr
  betabar_new = betabar
  index = 1:Dstar
  for (d in 1:Dstar){
    index_notd = index[-d]
    xnotd = betabar_new[index_notd]
    
    betabar_d = mean_betabar_update[d]
    betabar_notd = mean_betabar_update[index_notd]
    
    cov_dd = Sigmabetabar_update[d,d]
    cov_notdnotd = Sigmabetabar_update[index_notd, index_notd]
    cov_notdd = Sigmabetabar_update[cbind(rep(d, length(index_notd)), index_notd)]
    
    A = crossprod(as.matrix(cov_notdd), chol2inv(chol(cov_notdnotd)))
    con_meand = betabar_d + A %*% as.matrix(xnotd - betabar_notd) 
    con_sigmad = cov_dd - A %*% as.matrix(cov_notdd)
    
    u = runif(1, 0, 1)
    lb = parHP$betabar_par[d, 1]
    ub = parHP$betabar_par[d, 2]
    Fz = pnorm(lb, con_meand, con_sigmad) + u*(pnorm(ub, con_meand, con_sigmad) - pnorm(lb, con_meand, con_sigmad))
    
    if(Fz != 1 & Fz != 0)
    {
      betabar_new[d] = qnorm(Fz, con_meand, con_sigmad)
    } else {
      betabar_new[d] = parHP$betabar_par[d, Fz+1]
    }
  }
  
  #  Update parameters cond poster upsilonbeta
  shape_upsilonbeta = 1 + nubeta*Dstar/2
  rate_upsilonbeta = nubeta*(tr(chol2inv(chol(B_new))))/2
  
  #  Sample upsilonbeta from inverse cdf method on truncated distr
  u  = runif(1, 0, 1)
  lb = upsilonbeta_par[1]
  ub = upsilonbeta_par[2]
  Fz = pgamma(lb, shape = shape_upsilonbeta, rate = rate_upsilonbeta) + u*(pgamma(ub, shape = shape_upsilonbeta, rate = rate_upsilonbeta) - pgamma(lb,  shape = shape_upsilonbeta, rate = rate_upsilonbeta))
  if(Fz != 1 & Fz != 0)
  {
    upsilonbeta_new = qgamma(Fz,  shape = shape_upsilonbeta, rate = rate_upsilonbeta)
    
  } else {
    upsilonbeta_new = parHP$upsilonbeta_par[Fz+1]
  }
  
  return(list(mu_mu_new = mu_mu_new, Sigma_mu_new = Sigma_mu_new, upsilon_mu_new = upsilon_mu_new, upsilon_Sigma_new = upsilon_Sigma_new, a0_new = a0_new, betabar_new = betabar_new, B_new = B_new, upsilonbeta_new = upsilonbeta_new))
}
