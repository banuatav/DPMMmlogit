dpmm_mnl_gibbs = function(y, y_p, X_con, X_cat, X_p, Dummy, Dummy_p, PriorPar, MCMCpar, I = NULL){
  ################################## General ######################################
  library(mvtnorm)
  library(MCMCpack)
  library(psych) 
  library(dummies)
  if(is.null(I)){
    
    start_sample_number = MCMCpar$start_sample_number
    start_iter_number = MCMCpar$start_iter_number
    number_samples = MCMCpar$number_samples
    samplesize     = MCMCpar$samplesize
    m              = MCMCpar$m
    maxClus        = MCMCpar$maxClus
    AcceptRate     = rep(0, maxClus)
    ClusterOpen    = rep(0, maxClus)
    s              = MCMCpar$s
    thin           = MCMCpar$thin
    burnin         = MCMCpar$burnin
    
    J     = length(unique(y))
    n     = length(y)
    n_p   = length(y_p)
    Dcon  = dim(X_con)[2]
    Dcat  = dim(X_cat)[2]
    M_d   = apply(X_cat, 2, max)
    maxMd = max(M_d)
    X     =  cbind(X_con, X_cat)
    Dstar = (J-1)*(1+Dcon+sum(M_d-1))
    
    ################################## Create Space Parameters #####################################
    # Space parameters
    muStar      = matrix(NA, nrow = Dcon, ncol  = maxClus) # ~N(mu_0, sigma_0)
    SigmaStar   = matrix(NA, nrow = Dcon*Dcon, ncol  = maxClus) # ~SIX2(nu,tau2)
    piStar      = matrix(NA, nrow = Dcat*maxMd, ncol = maxClus) # ~Dir(a_0/Md)
    betaStar    = matrix(NA, nrow = Dstar, ncol  = maxClus) # ~N(betabar, B)
    alpha       = NA # ~Ga(a,b)
    
    ################################## Initialize Parameters #####################################
    # Hyperparameters
    mu_mu = rowMeans(PriorPar$mu_mu_par)
    upsilon_mu = mean(PriorPar$upsilon_mu_par)
    nu_mu = PriorPar$nu_mu
    Sigma_mu  = riwish_rcpp(nu_mu, nu_mu*upsilon_mu*diag(rep(1,Dcon))) # IW(nu_mu, nu_mu*upsilon_mu*I)
    
    nu_Sigma = PriorPar$nu_Sigma
    upsilon_Sigma = mean(PriorPar$upsilon_Sigma_par)  
    
    a0              = rowMeans(PriorPar$a0_par)
    
    upsilonbeta    = mean(PriorPar$upsilonbeta_par)
    nubeta         = PriorPar$nubeta
    betabar        = rowMeans(PriorPar$betabar_par)
    B              = riwish_rcpp(nubeta, nubeta*upsilonbeta*diag(rep(1,Dstar))) # IW(nubeta, nubeta*upsilonbeta*I)
    
    # Cluster assignment, assigmnent initialized as all belong to cluster 1
    if(is.null(PriorPar$z)){ z = rep(1, n) } else {z = PriorPar$z}
    N_k = rep(NA, maxClus)
    K   = length(unique(z))
    N_k[1:K] = as.vector(table(z))
    
    # Parameters 
    for(k in 1:K){
      SigmaStar[, k]    = c(riwish_rcpp(nu_Sigma, nu_Sigma*upsilon_Sigma*diag(1,Dcon))) # ~IW(nu_Sigma, nu_Sigma * upsilon_Sigma*I)
      muStar[, k]       = rmvnormrcpp(1, mu_mu, Sigma_mu) # ~N(mu_mu, Sigma_mu)  
      betaStar[,k]      = rmvnormrcpp(1, betabar, B) # ~N(betabar, B)
      for (d in 1:Dcat){ piStar[((d-1)*maxMd+1):((d-1)*maxMd+M_d[d]), k] = rdirichlet(n = 1, alpha = (rep(a0[d], M_d[d])/M_d[d])) } # ~Dir(a0/M_d) 
    }
    
    alpha           = 0.0001
    
    thetaStar = list(betaStar = betaStar, muStar=muStar, SigmaStar=SigmaStar, piStar=piStar, alpha = alpha)
    hyperPar = list(mu_mu = mu_mu, upsilon_mu = upsilon_mu, nu_mu = nu_mu, Sigma_mu = Sigma_mu, nu_Sigma = nu_Sigma, upsilon_Sigma = upsilon_Sigma, a0 = a0,  betabar = betabar, B = B, upsilonbeta = upsilonbeta, nubeta = nubeta)
    parHP = list( mu_mu_par = PriorPar$mu_mu_par, upsilon_mu_par = PriorPar$upsilon_mu_par, upsilon_Sigma_par = PriorPar$upsilon_Sigma_par, a0_par = PriorPar$a0_par, betabar_par = PriorPar$betabar_par, upsilonbeta_par = PriorPar$upsilonbeta_par, alpha_par =  PriorPar$alpha_par)
    testSample = list(y_pred = y_p, X_pred = X_p, Dummy_pred = Dummy_p, n_pred = n_p)
    DAInfo = list(y = y, Dy = dummy(y),  X = X, Dummy = Dummy, n = n,  Dcon = Dcon, Dcat = Dcat, Dstar = Dstar, J = J, M_d = M_d, maxMd = maxMd, m = m, AcceptRate = AcceptRate,   ClusterOpen =   ClusterOpen, s=s)
    clusterInfo = list(z = z, K = K, N_k = N_k)
    remove(betaStar, muStar, SigmaStar, piStar, alpha, a0, betabar, B, upsilonbeta, nubeta, mu_mu, upsilon_mu, nu_mu, Sigma_mu, nu_Sigma, upsilon_Sigma)
    remove(y_p, X_p, Dummy_p, n_p, y, X, Dummy, n, Dcon, Dcat, Dstar, J, M_d, maxMd, m, AcceptRate, ClusterOpen, s, z, K, N_k)
    remove(PriorPar, MCMCpar)
    sample_iteration_index = start_sample_number:number_samples
    samplesize_iteration_index = start_iter_number:samplesize
  }else{
    
    thetaStar   = I$thetaStar
    clusterInfo = I$clusterInfo
    hyperPar    = I$hyperPar
    
    m               = MCMCpar$m
    number_samples  = MCMCpar$number_samples
    samplesize      = MCMCpar$samplesize
    maxClus         = MCMCpar$maxClus
    
    start_sample_number = MCMCpar$start_sample_number
    start_iter_number   = MCMCpar$start_iter_number
    
    sample_iteration_index     = start_sample_number:number_samples
    samplesize_iteration_index = start_iter_number:samplesize
    
    AcceptRate     = rep(0, maxClus)
    ClusterOpen    = rep(0, maxClus)
    s              = MCMCpar$s
    thin           = MCMCpar$thin
    burnin         = MCMCpar$burnin
    
    J     = length(unique(y))
    n     = length(y)
    n_p   = length(y_p)
    Dcon  = dim(X_con)[2]
    Dcat  = dim(X_cat)[2]
    M_d   = apply(X_cat, 2, max)
    maxMd = max(M_d)
    X     =  cbind(X_con, X_cat)
    Dstar = (J-1)*(1+Dcon+sum(M_d-1))
    
    testSample = list(y_pred = y_p, X_pred = X_p, Dummy_pred = Dummy_p, n_pred = n_p)
    DAInfo = list(y = y, Dy = dummy(y),  X = X, Dummy = Dummy, n = n,  Dcon = Dcon, Dcat = Dcat, Dstar = Dstar, J = J, M_d = M_d, maxMd = maxMd, m = m, AcceptRate = AcceptRate,   ClusterOpen =   ClusterOpen, s=s)
    parHP = list( mu_mu_par = PriorPar$mu_mu_par, upsilon_mu_par = PriorPar$upsilon_mu_par, upsilon_Sigma_par = PriorPar$upsilon_Sigma_par, a0_par = PriorPar$a0_par, betabar_par = PriorPar$betabar_par, upsilonbeta_par = PriorPar$upsilonbeta_par, alpha_par =  PriorPar$alpha_par)
    
    remove(y_p, X_p, Dummy_p, n_p)
    remove(y, X, Dummy, n, Dcon, Dcat, Dstar, J, M_d, maxMd, m, AcceptRate,   ClusterOpen, s)
    remove(I, MCMCpar, d)
    
  }
  
  ################################## MCMC Iterations ######################################
  
  n           = DAInfo$n
  Dcon        = DAInfo$Dcon
  Dcat        = DAInfo$Dcat
  Dstar       = DAInfo$Dstar
  J           = DAInfo$J 
  maxMd       = DAInfo$maxMd 
  M_d         = DAInfo$M_d
  K           = clusterInfo$K
  
  predictionprob = matrix(0, nrow = testSample$n_pred, ncol = 2 * J + 1)
  llKeep = matrix(NA, ceiling(samplesize/thin), 4)
  alphaKeep = rep(NA, ceiling(samplesize/thin))
  KKeep = rep(NA, ceiling(samplesize/thin))
  zKeep = matrix(NA, ceiling(samplesize/thin), n)
  betaKeep  = array(NA, c(ceiling(samplesize/thin), Dstar, maxClus))
  SigmaKeep = array(NA, c(ceiling(samplesize/thin), Dcon*Dcon, maxClus))
  muKeep    = array(NA, c(ceiling(samplesize/thin), Dcon, maxClus))
  piKeep    = array(NA, c(ceiling(samplesize/thin), Dcat*maxMd, maxClus))
  mu_muKeep  = matrix(NA, ceiling(samplesize/thin), Dcon)
  ups_muKeep = rep(NA, ceiling(samplesize/thin))
  Sig_muKeep = matrix(NA, ceiling(samplesize/thin), Dcon*Dcon)
  ups_SigKeep = rep(NA, ceiling(samplesize/thin))
  a0Keep = matrix(NA, ceiling(samplesize/thin), Dcat)
  betabarKeep = matrix(NA, ceiling(samplesize/thin), Dstar)
  ups_betaKeep = rep(NA, ceiling(samplesize/thin))
  BKeep =  matrix(NA, ceiling(samplesize/thin), Dstar*Dstar)
  
  for (sample_nr in sample_iteration_index){
    
    DAInfo$AcceptRate     = rep(0, maxClus)
    DAInfo$ClusterOpen    = rep(0, maxClus)
    
    for(j in samplesize_iteration_index){
      # Sample Clusterassignments
      if(j%%10 -1 == 0){
        cat('Sample =', sample_nr, '. Iteration =', j, '. ', '#Clusters =', K, '. ', fill = TRUE)
      }
      
      result = sampleCluster(thetaStar = thetaStar, clusterInfo = clusterInfo , DAInfo = DAInfo, hyperPar = hyperPar)
      clusterInfo = result$clusterInfo
      thetaStar   = result$thetaStar
      K           = clusterInfo$K
      
      resultTheta = sampleTheta(thetaStar = thetaStar, clusterInfo = clusterInfo , DAInfo = DAInfo, hyperPar = hyperPar)
      thetaStar   = resultTheta$thetaStar
      DAInfo$AcceptRate = resultTheta$AcceptRate
      DAInfo$ClusterOpen = resultTheta$ClusterOpen
      
      if(j%%50 - 1 == 0){
        cat("Acceptrate = ", (DAInfo$AcceptRate[1:K]/DAInfo$ClusterOpen[1:K]), fill = TRUE)
        cat("Clustersizes = ", clusterInfo$N_k[1:K], fill = TRUE)
      }
      
      # 3. Sample alpha
      resultAlpha             = sampleAlpha(alpha_current = thetaStar$alpha, n = n, K = K, a = parHP$alpha_par[1], b = parHP$alpha_par[2])
      thetaStar$alpha         = resultAlpha$alpha_new
      palphaKeep              = resultAlpha$p_alpha
      
      #4. Sample hyperparameters
      resultHyp              = sampleHyp(parHP = parHP, thetaStar = thetaStar, hyperPar = hyperPar, K = K , M_d = M_d, maxMd = maxMd, Dcon = Dcon,  Dcat = Dcat,  Dstar = Dstar)
      hyperPar$mu_mu         = resultHyp$mu_mu_new
      hyperPar$upsilon_mu    = resultHyp$upsilon_mu_new
      hyperPar$Sigma_mu      = resultHyp$Sigma_mu_new
      hyperPar$upsilon_Sigma = resultHyp$upsilon_Sigma_new
      
      hyperPar$a0            = resultHyp$a0_new
      
      hyperPar$betabar       = resultHyp$betabar_new
      hyperPar$B             = resultHyp$B_new
      hyperPar$upsilonbeta   = resultHyp$upsilonbeta_new
      
      
      if(j%%thin == 0 ){
        llKeep[floor(j/thin),]          = CCmeasure(DAInfo, thetaStar, clusterInfo, hyperPar, parHP$alpha_par)
        alphaKeep[floor(j/thin)]        = thetaStar$alpha
        KKeep[floor(j/thin)]            = K
        zKeep[floor(j/thin),]           = clusterInfo$z
        
        betaKeep[floor(j/thin), , 1:K]  = thetaStar$betaStar[,1:K]
        SigmaKeep[floor(j/thin), , 1:K] = thetaStar$SigmaStar[,1:K]
        muKeep[floor(j/thin), , 1:K]    = thetaStar$muStar[,1:K]
        piKeep[floor(j/thin), , 1:K]    = thetaStar$piStar[,1:K]
        
        mu_muKeep[floor(j/thin), ]    =  hyperPar$mu_mu
        ups_muKeep[floor(j/thin)]     =  hyperPar$upsilon_mu
        Sig_muKeep[floor(j/thin), ]   =  c(hyperPar$Sigma_mu)
        
        ups_SigKeep[floor(j/thin)]    = hyperPar$upsilon_Sigma
        
        a0Keep[floor(j/thin), ]       = hyperPar$a0
        
        betabarKeep[floor(j/thin), ]  = hyperPar$betabar
        ups_betaKeep[floor(j/thin)]   = hyperPar$upsilonbeta
        BKeep[floor(j/thin), ]        =  c(hyperPar$B)
        
        if(sample_nr>burnin){
          predictionprob    = predictionprob + prediction(testSample = testSample,  thetaStar = thetaStar, K_t =  K, N_k = clusterInfo$N_k, n = n, J = DAInfo$J, M_d = M_d, maxMd = maxMd, Dcon = DAInfo$Dcon, Dcat = DAInfo$Dcat, Dstar = DAInfo$Dstar, hyperPar,  m = precisionmeasure)
        }
      }
      
    }
    
    name_file_I = paste("I", sample_nr,"_" , j ,".RData", sep = "")
    I_s = list(sample_nr = sample_nr, iter = j, thetaStar = thetaStar, clusterInfo = clusterInfo , hyperPar = hyperPar)
    save(I_s, file = name_file_I)
    
    CCdiag = list(ll = llKeep, alpha = alphaKeep)
    save(CCdiag, file = paste0("CCdiag", sample_nr,".RData"))
    
    SaveThetaStar = list(beta = betaKeep, mu = muKeep, Sigma = SigmaKeep, pi = piKeep)
    save(SaveThetaStar, file = paste0("ThetaStar", sample_nr,".RData"))
    
    SaveClusterInfo = list(K = KKeep, z = zKeep)
    save(SaveClusterInfo, file = paste0("ClusterInfo", sample_nr,".RData"))
    
    SaveHP = list(mu_mu = mu_muKeep, ups_mu = ups_muKeep, Sig_mu = Sig_muKeep, ups_Sig = ups_SigKeep, a0 = a0Keep, betabar = betabarKeep, ups_beta = ups_betaKeep, B = BKeep)
    save(SaveHP, file = paste0("HP", sample_nr,".RData"))
    
    if(sample_nr>burnin){
      save(predictionprob, file = paste0("PredictionProb", sample_nr,".RData"))
    }
    samplesize_iteration_index = 1:samplesize
    
  }
  return(number_samples) 
}
