sampleAlpha = function(alpha_current, n, K, a, b){
  ##### Samples parameter alpha in a Dirichlet process G~DP(alpha, G_0)
  ##### under prior specification alpha ~ Ga(a, b)
  ##### auxiliary variable eta is samples along to make this Gibbs step possible
  #### inputs: alpha_current = current value of alpha, n = number of observations, K = current number of clusters
  
  # First sample auxiliary variable eta from ~Beta(alpha + 1, n)
  eta = rbeta(n = 1, shape1 = alpha_current + 1, shape2 = n)
  
  # Then sample alpha given eta. The distrbution of alpha given eta is a mixture of two Gammas with
  # mixing constant c_eta. 
  c_eta = (a + K - 1)/(a + K - 1 + n*(b - log(eta)))
  u = runif(n = 1, min = 0, max =1)
  if(u<c_eta){
    alpha_new  = rgamma(1, shape = a + K, rate  = b - log(eta)) 
    p_alpha    = pgamma(alpha_new, shape = a + K, rate =  b - log(eta))
  }else{
    alpha_new  = rgamma(1, shape = a + K - 1, rate = b - log(eta))
    p_alpha    = pgamma(alpha_new, shape = a + K -1, rate =  b - log(eta))
  }
  
  return(list(alpha_new = alpha_new, p_alpha = p_alpha))
}