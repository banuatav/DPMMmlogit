#include <RcppArmadillo.h>

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
// [[Rcpp::export]]
arma::vec log_dmvnorm_arma_mult(arma::mat x, arma::mat mean, arma::mat sigma) { 
  int n = x.n_rows;
  int dim = x.n_cols;
  arma::mat x_cen(n, dim);
  arma::mat cov(dim,dim);
  arma::vec logdet(n);
  arma::mat x_stan(n, dim);
  x_cen = x - mean;
  
  for (int i=0; i < n; i++) {
    
    for (int c = 0; c < dim; c++) {
      for (int r = 0; r < dim; r++) {
        cov(r,c) = sigma(i,c*dim + r);
      }
    }
    x_stan.row(i) =  x_cen.row(i) * cov.i();
    logdet(i) = sum(arma::log(arma::eig_sym(cov)));
  }
  
  arma::vec logretval = -((dim*log2pi + logdet + sum(x_stan % x_cen, 1))/2  ) ;
  

  return(logretval);

}
// [[Rcpp::export]]
double log_dmvnorm_arma_sing(arma::vec x, arma::vec mean, arma::mat sigma) { 

  int dim = sigma.n_cols;
  arma::vec x_cen;
  x_cen = x - mean;
  arma::mat xSigx = trans(x_cen) * sigma.i() * x_cen;
  double sumsq;
  sumsq = as_scalar(xSigx);
  double logdens;
  logdens = -((dim*log2pi + sum(arma::log(arma::eig_sym(sigma))) + sumsq)/2  ) ;

  return(logdens);
  
}
// [[Rcpp::export]]
arma::mat invHess(arma::mat B, arma::vec y_k, arma::mat R_k, arma::mat betaMat, int J){
  arma::mat Bi = B.i();
  arma::mat H_k;
  
  int n = R_k.n_rows;
  int dim = R_k.n_cols;
  int Dstar = dim*(J-1);
  
  arma::mat Xbeta;
  Xbeta = R_k * betaMat;
  arma::vec maxXbeta = Xbeta.col(0);

  for (int r = 0; r<n; r++){
    for (int c = 1; c<J; c++){
      if(Xbeta(r,c)>maxXbeta(r)){
        maxXbeta(r) = Xbeta(r,c);
      }
    }
  }


  arma::rowvec iota(J);
  iota.fill(1);

  Xbeta = exp(Xbeta - maxXbeta * iota);
  arma::vec sumXbeta = sum(Xbeta,1);
  arma::mat denom = sumXbeta * iota;
  arma::mat Prob = Xbeta/denom;
  arma::mat Hess(Dstar, Dstar);


  Hess.fill(0);
  arma::mat Xt;
  arma::mat Diag;
  arma::vec ones(J-1);
  ones.fill(1);
  Diag = diagmat(ones);
  arma::mat zeros(1,Dstar);
  zeros.fill(0);
  
  
  arma::vec p(J);
  arma::mat A(J,J);
  arma::mat X_t;

  for (int i = 0; i<n; i++) {
    p = trans(Prob.row(i));
    A = diagmat(p) - p * trans(p);
    Xt =  join_cols(kron(R_k.row(i), Diag), zeros);
    Hess = Hess +  trans(Xt) * A * Xt;
    }

  H_k = Hess + Bi;
  
  return(H_k.i());
}
// [[Rcpp::export]]
arma::mat rmvnormrcpp(int n_, arma::vec mu_, arma::mat sigma_){
  int n = n_;
  arma::vec mu = mu_;
  arma::mat sigma = sigma_;
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return (repmat(mu, 1, n).t() + Y *chol(sigma));
}
// [[Rcpp::export]]
arma::mat rmvnormrcpp_chol(int n, arma::vec mu, arma::mat col_sigma){
  int ncols = col_sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return (repmat(mu, 1, n).t() + Y *col_sigma);
}
// [[Rcpp::export]]
arma::mat chol_rcpp(arma::mat matrix){
  return (chol(matrix));
}
// [[Rcpp::export]]
arma::mat riwish_rcpp(double nu, arma::mat upsilon){
  int dim = upsilon.n_cols;
  arma::mat upsilon_wis = upsilon.i();
  arma::vec mu(dim);
  mu.fill(0);
  arma::mat Y = arma::randn(nu, dim);
  arma::mat X = (repmat(mu, 1, nu).t() + Y *chol(upsilon_wis));
  return (X.t() * X).i();
}
// [[Rcpp::export]]
arma::mat riwish_rcpp_chol(double nu, arma::mat chol_inverse_upsilon){
  int dim = chol_inverse_upsilon.n_cols;
  arma::vec mu(dim);
  mu.fill(0);
  arma::mat Y = arma::randn(nu, dim);
  arma::mat X = (repmat(mu, 1, nu).t() + Y *chol_inverse_upsilon);
  return (X.t() * X).i();
}
// [[Rcpp::export]]
double logdiwish_rcpp(double nu, double upsilon, arma::mat Draw){
  int dim = Draw.n_cols;
  arma::mat DrawInv = Draw.i();
  double detDraw = det(Draw);
  return (nu*dim/2 * log(nu*upsilon)-(nu+dim+1)/2*log(detDraw)-nu*upsilon*trace(DrawInv)/2);
}
// [[Rcpp::export]]
arma::mat prob_mlogit(arma::mat beta, arma::mat R_p){
  int n_p = R_p.n_rows;
  int J = beta.n_cols;
  arma::mat Xbeta = arma::zeros(n_p,J+1);
  Xbeta.submat(0, 0, n_p-1, J-1) = R_p * beta;
  arma::vec maxXbeta = max(Xbeta, 1);
  arma::mat relXbeta     = Xbeta - repmat(maxXbeta,1,J+1);
  arma::vec denom        = log(sum(exp(relXbeta),1));
  return(relXbeta  - repmat(denom,1,J+1) );
}
