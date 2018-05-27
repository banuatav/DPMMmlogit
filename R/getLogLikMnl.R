getLogLikMnl = function (beta, Dy, X, J) {
  n = nrow(Dy)
  Dstar = ncol(X) * (J-1)
  Xbeta = cbind(X %*% matrix(beta, nrow = Dstar/(J-1), ncol = J-1, byrow = TRUE), rep(0, n))
  Xbeta = Xbeta - apply(Xbeta, 1, max)
  denom = rowSums(exp(Xbeta))
  lProb = Xbeta- log(denom)
  ll = sum(lProb[Dy==1])
  return(ll)
}