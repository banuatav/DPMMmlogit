logGammaApprox = function(z){
  return((log(2*pi)-log(z))/2 + z*(log(z+1/(12*z -1/10/z))-1))
}
CCmeasure = function(DAInfo, thetaStar, clusterInfo, hyperPar, alpha_par){
  y       = DAInfo$y
  X       = DAInfo$X
  Dummy   = DAInfo$Dummy
  Dy      = DAInfo$Dy
  J       = DAInfo$J
  M_d     = DAInfo$M_d
  maxMd   = DAInfo$maxMd
  Dcon    = DAInfo$Dcon
  Dcat    = DAInfo$Dcat
  Dstar   = DAInfo$Dstar
  n       = DAInfo$n
  
  nu_beta = hyperPar$nubeta
  nu_mu = hyperPar$nu_mu
  nu_Sigma = hyperPar$nu_Sigma
  
  z_j     = clusterInfo$z
  alpha_j = thetaStar$alpha
  K_j     = clusterInfo$K
  beta_j  = thetaStar$beta
  Sigma_j = thetaStar$Sigma
  mu_j    = thetaStar$mu
  pi_j    = thetaStar$pi
  
  mu_mu_j      = hyperPar$mu_mu 
  upsilon_mu_j = hyperPar$upsilon_mu 
  Sigma_mu_j    = hyperPar$Sigma_mu
  
  upsilon_Sigma_j = hyperPar$upsilon_Sigma
  
  a0_j = hyperPar$a0
  
  betabar_j      = hyperPar$betabar
  upsilon_beta_j = hyperPar$upsilonbeta
  B_j =  hyperPar$B
  ll = rep(NA,4)
  ll_c_j  = 0
  ll_thetaStar_j = 0
  for(k in 1:K_j){
    ## Select Data
    y_k  = y[z_j==k]
    n_k   = length(y_k)
    if ( n_k == 1 ) {
      X_k         = t(as.matrix(X[z_j==k, ]))
      X_con_k     = t(as.matrix(X_k[,1:Dcon]))
      X_cat_k     = t(as.matrix(X_k[,-(1:Dcon)]))
      Dummy_k     = t(as.matrix(Dummy[z_j == k, ]))
      Dy_k        = t(as.matrix(Dy[z_j == k, ]))
    } else {
      X_k  = X[z_j==k, ] 
      X_con_k   = as.matrix(X_k[,1:Dcon])
      X_cat_k     = as.matrix(X_k[,-(1:Dcon)])
      Dummy_k     = Dummy[z_j == k, ]
      Dy_k        = Dy[z_j == k, ]
    }
    R_k           = cbind(rep(1, n_k), X_con_k, Dummy_k)
    
    mu_j_k        = as.matrix(mu_j[,k])
    Sigma_j_k     = matrix(Sigma_j[,k], nrow = Dcon)
    pi_j_k        = matrix(pi_j[,k], ncol = Dcat)
    beta_j_k      = as.matrix(beta_j[,k])
    beta_matrix   = matrix(beta_j_k, ncol = (J-1), nrow = length(beta_j_k)/(J-1), byrow = TRUE)
    
    Xbeta      = cbind(R_k %*% beta_matrix, rep(0, n_k))
    maxXbeta   = apply(Xbeta, 1, max)
    relXbeta   = Xbeta - maxXbeta  #MASimpl matrix - vector
    denom      = log(rowSums(exp(relXbeta)))
    ll_y_x_ll  = sum(relXbeta[cbind(1:n_k, y_k)] - denom)  #MASimpl matrix - vector
    ll_xcon_ll = sum(dmvnorm(x = X_con_k, mean = mu_j_k, sigma = Sigma_j_k, log = TRUE))
    ll_xcat_ll = 0
    for (d in 1:Dcat){
      ll_xcat_ll = ll_xcat_ll + sum(as.matrix(log(pi_j_k[cbind(X_cat_k[,d], rep(d, length(X_cat_k[,d])))])))
    }
    
    ll_xy = ll_y_x_ll + ll_xcon_ll + ll_xcat_ll
    
    ll_pi = 0 
    for (d in 1:Dcat){
      a_d = a0_j[d]/M_d[d]
      #ll_pi = ll_pi + log(ddirichlet(pi_j_k[1:M_d[d],d], a_d))
      if(a0_j[d]<170){
        ll_pi = ll_pi + sum(log(pi_j_k[1:M_d[d],d]))*(a_d-1) + log(gamma(a0_j[d])) - M_d[d] *log(gamma(a_d[1]))
      } else if(a0_j[d]>170 & a_d<170){
        ll_pi = ll_pi + sum(log(pi_j_k[1:M_d[d],d]))*(a_d-1) + logGammaApprox(a0_j[d]) - M_d[d] *log(gamma(a_d[1]))
      }else{
        ll_pi = ll_pi + sum(log(pi_j_k[1:M_d[d],d]))*(a_d-1) + logGammaApprox(a0_j[d]) - M_d[d] *logGammaApprox(a_d[1])
      }
    }
    ll_thetaStar_j = ll_thetaStar_j + ll_pi + dmvnorm(c(beta_j_k), betabar_j, B_j, log = TRUE) + dmvnorm(c(mu_j_k), mu_mu_j, Sigma_mu_j, log = TRUE) + log(diwish(Sigma_j_k, nu_Sigma, nu_Sigma*upsilon_Sigma_j*diag(Dcon)))
    ll_c_j         =  ll_c_j  + ll_xy
    
  }
  ll[1] =  ll_c_j 
  ll[2] =  ll_thetaStar_j
  ll[3] =  logdiwish_rcpp(nu_beta, upsilon_beta_j, B_j)
  ll[4] =  logdiwish_rcpp(nu_mu, upsilon_mu_j, Sigma_mu_j)
  return(ll)
}