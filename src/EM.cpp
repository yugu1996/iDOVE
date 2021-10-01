#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <ctime>
#include "BS.h"

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

/* EM algorithm when covariates are present in the model

  w	  inout	{n x nTimes}
  beta	  out	{nBeta}
  gamma   out	{nGam}
  lambda  out   {nTimes} 
  L       in    {n} left interval times
  R       in    {n} right interval times
  Rstar   in    {n} right interval times modified to remove Inf
  t       in    {nTimes} time points
  tind    in    {n x 3} indices for time nearest entry, left, and right times
  X       in    {n x nBeta} covariates of model
  S       in    {n} Vaccination time
  knots   in    {nKnots} knots for B-spline
  constantVE in {1} indicator if slope of last spline is zero
  threshold in  {1} value below which algorithm is deemed to have converged
  maxit   in    {1} maximum number of iterations to consider
*/

// [[Rcpp::export]]
void EM(arma::mat& w, arma::vec& beta, arma::vec& gamma, arma::vec& lambda, 
        const arma::ivec& L, const arma::ivec& R, const arma::ivec& Rstar, 
        const arma::ivec& t, const arma::imat& tind, const arma::mat& X, 
        const arma::ivec& S, const arma::ivec& knots, bool constantVE, 
        double threshold=10^(-4), int maxit=5000) {

  int n = L.n_elem, m = t.n_elem, p1 = beta.n_elem, p2 = gamma.n_elem;
  double s, exps, exbeta, texbeta, dif = 1.0;
  arma::vec beta0, gamma0, lambda0, temp, U, Uvec, bvec; 
  arma::rowvec Sw, S0, tmpe;
  arma::colvec XZ;
  arma::mat I, I_inv, S1, Swx;
  arma::cube S2;
  bool flag;

  bvec.zeros(m);

  int it = 0;

  while (it<=maxit && dif>threshold) {
    it++;

    R_CheckUserInterrupt();

    beta0 = beta; 
    gamma0 = gamma; 
    lambda0 = lambda;

    // begin: E-step
    for (int i=0; i<n; i++) {

      if (tind(i,1)+1 > tind(i,2)) continue;

      exbeta = exp(arma::as_scalar(X.row(i)*beta));
      s = 0.0; 
      for (int k=tind(i,1)+1; k<=tind(i,2); k++) {
        if (S(i)<t(k)) {
          bvec(k) = lambda(k)*exp(arma::as_scalar(BS(t(k)-S(i), knots, constantVE)*gamma));
        } else {
          bvec(k) = lambda(k);
        }
        s += bvec(k);
      }
      s *= exbeta;

      if (s>1e-16) {
        exps = exbeta / (1.0 - exp(-s));
        for (int k=tind(i,1)+1; k<=tind(i,2); k++) {
          w(i,k) = bvec(k)*exps;
        }
      }
    }
    // end: E-step
    // begin: M-step

    //// begin: update theta

    U = arma::zeros<arma::vec>(p1+p2);
    I = arma::zeros<arma::mat>(p1+p2, p1+p2);
    S0.zeros(m); 
    S1.zeros(p1+p2, m); 
    S2.zeros(p1+p2, p1+p2, m); 
    Sw.zeros(m); 
    Swx.zeros(p1+p2, m);
    int i = n-1;
    for (int k=m-1; k>=0; k--) {
      while (i>=0 && Rstar(i)>=t(k)) {
        texbeta = exp(arma::as_scalar(X.row(i)*beta));
        for(int l=k; l>=tind(i,0); --l) {

          if (S(i)<t(l)) {
            tmpe = BS(t(l)-S(i), knots, constantVE);
            exbeta = texbeta*exp(arma::as_scalar(tmpe*gamma));
          } else {
            tmpe.zeros(p2);
            exbeta = texbeta;
          }
          XZ = trans(join_horiz(X.row(i), tmpe));

          S0(l) += exbeta;
          S1.col(l) += exbeta*XZ;
          S2.slice(l) += exbeta*XZ*trans(XZ);
          Sw(l) += w(i,l);
          Swx.col(l) += XZ*w(i,l);
        }
        i--;
      }
      Uvec = S1.col(k)/S0(k);
      U += Swx.col(k) - Uvec*Sw(k);
      I += Sw(k)*(S2.slice(k)/S0(k)-Uvec*Uvec.t());
    }
    flag = pinv(I_inv, I);

    if (!flag) {
      stop("EM terminated due to singular information matrix");
    }

    temp = I_inv*U;
    beta += temp.head(p1); 
    gamma += temp.tail(p2);
    //// end: update theta

    //// begin: update lambda
    i = n-1;
    S0.zeros(m); 
    Sw.zeros(m);
    for (int k=m-1; k>=0; k--) {
      while (i>=0 && Rstar(i)>=t(k)) {
        exbeta = exp(arma::as_scalar(X.row(i)*beta));
        for (int l=k; l>=tind(i,0); --l) {
          if (S(i)<t(l)) {
            S0(l) += exbeta*exp(arma::as_scalar(BS(t(l)-S(i), knots, constantVE)*gamma));
          } else {
            S0(l) += exbeta;
          }
          Sw(l) += w(i,l);
        }
        i--;
      }
      lambda(k) = Sw(k)/S0(k);
    }
    //// end: update lambda

    // end: M-step
    dif = norm(beta-beta0, 2) + norm(gamma-gamma0, 2) + 
          norm(lambda-lambda0, "inf");
    if (it % 100 == 0) {
      Rprintf("Iteration %u : difference = %f \n", it, dif);
    }

  }

  if (dif>threshold) {
    warning("EM algorithm did not converge");
  } else {
    Rprintf("EM algorithm converged after %u iterations \n", it);
  }

}

/* log likelihood when covariates are present in the model

  beta	  out	{nBeta}
  gamma   out	{nGam}
  lambda  out   {nTimes} 
  L       in    {n} left interval times
  R       in    {n} right interval times
  t       in    {nTimes} time points
  tind    in    {n x 3} indices for time nearest entry, left, and right times
  X       in    {n x nBeta} covariates of model
  S       in    {n} Vaccination time
  knots   in    {nKnots} knots for B-spline
  constantVE in {1} indicator if slope of last spline is zero
*/

// [[Rcpp::export]]
arma::vec LL(const arma::vec& beta, const arma::vec& gamma, 
             const arma::vec& lambda, const arma::ivec& L, const arma::ivec& R, 
             const arma::ivec& t,const arma::imat& tind, const arma::mat& X, 
             const arma::ivec& S, const arma::ivec& knots, bool constantVE) {
  int n = L.n_elem;
  double SLi, SUi, exbeta;
  arma::vec ll(n);

  for (int i=0; i<n; ++i) {
    exbeta = exp(arma::as_scalar(X.row(i)*beta));

    SLi = 0.0;
    for (int k=tind(i,0); k<=tind(i,1); k++) {
      if (S(i)<t(k)) {
        SLi += lambda(k) * exbeta *
               exp(arma::as_scalar(BS(t(k)-S(i), knots, constantVE)*gamma));
      } else {
        SLi += lambda(k) * exbeta;
      }
    }

    if (L(i)==R(i)) {
      if (S(i)<t(tind(i,1))) {
        ll(i) = log(lambda(tind(i,1)) * exbeta * exp(-SLi) *
                    exp(arma::as_scalar(BS(t(tind(i,1))-S(i), knots, constantVE)*gamma)));
      } else {
        ll(i) = log(lambda(tind(i,1)) * exbeta * exp(-SLi));
      }
    } else if (R(i)>=0) {
      SUi = 0.0;
      for (int k=tind(i,1)+1; k<=tind(i,2); k++) {
        if (S(i)<t(k)) {
          SUi += lambda(k)*exbeta*exp(arma::as_scalar(BS(t(k)-S(i), knots, constantVE)*gamma));
        } else {
          SUi += lambda(k)*exbeta;
        }
      }
      SUi += SLi;
      ll(i) = log(exp(-SLi)-exp(-SUi));
    } else {
      ll(i) = -SLi;
    }
  }
  //return accu(ll);
  return ll;
}

/* EM algorithm for lambda

  w	  inout	{n x nTimes}
  beta	  in	{nBeta}
  gamma   in	{nGam}
  lambda_init  in   {nTimes} 
  L       in    {n} left interval times
  R       in    {n} right interval times
  Rstar   in    {n} right interval times modified to remove Inf
  t       in    {nTimes} time points
  tind    in    {n x 3} indices for time nearest entry, left, and right times
  X       in    {n x nBeta} covariates of model
  S       in    {n} Vaccination time
  knots   in    {nKnots} knots for B-spline
  constantVE in {1} indicator if slope of last spline is zero
  threshold in  {1} value below which algorithm is deemed to have converged
  maxit   in    {1} maximum number of iterations to consider
*/

// [[Rcpp::export]]
arma::vec PLL(arma::mat& w, const arma::vec& beta, const arma::vec& gamma, 
              const arma::vec& lambda_init, const arma::ivec& L, const arma::ivec& R, 
              const arma::ivec& Rstar, const arma::ivec& t, const arma::imat& tind, 
              const arma::mat& X, const arma::ivec& S, const arma::ivec& knots, 
              bool constantVE, double threshold=10^(-4), int maxit=5000) {
  int n = L.n_elem, m = t.n_elem, it;
  double s, exbeta, exps, dif = 1;
  arma::vec lambda0, bvec, lambda = lambda_init; 
  arma::rowvec Sw, S0, XZ;

  bvec.zeros(m);

  it = 0;
  while (it<maxit && dif>threshold) {

    R_CheckUserInterrupt();

    it++;

    lambda0 = lambda;

    // begin: E-step
    for (int i=0; i<n; ++i) {

      if (tind(i,2)<tind(i,1)+1) continue;

      exbeta = exp(arma::as_scalar(X.row(i)*beta));
      s = 0.0; 
      for (int k=tind(i,1)+1; k<=tind(i,2); k++) {
        if (S(i)<t(k)) {
          bvec(k) = lambda(k)*exp(arma::as_scalar(BS(t(k)-S(i), knots, constantVE)*gamma));
        } else {
          bvec(k) = lambda(k);
        }
        s += bvec(k);
      }
      s *= exbeta;

      if (s > 1e-16) {
        exps = exbeta / (1.0 - exp(-s));
        for (int k=tind(i,1)+1; k<=tind(i,2); k++) {
          w(i,k) = bvec(k)*exps;
        }
      }
    }
    // end: E-step

    // begin: M-step

    //// begin: update lambda
    int i = n-1;
    S0.zeros(m); 
    Sw.zeros(m);
    for (int k=m-1; k>=0; k--) {
      while (i>=0 && Rstar(i)>=t(k)) {
        for (int l=k; l>=tind(i,0); --l) {
          if (S(i)<t(l)) {
            S0(l) += exp(arma::as_scalar(X.row(i)*beta) +
                         arma::as_scalar(BS(t(l)-S(i), knots, constantVE)*gamma));
          } else {
            S0(l) += exp(arma::as_scalar(X.row(i)*beta));
          }
          Sw(l) += w(i,l);
        }
        i--;
      }
      lambda(k) = Sw(k)/S0(k);
    }
    //// end: update lambda

    // end: M-step

    dif = norm(lambda-lambda0, "inf");
  }

  if (dif>threshold) {
    warning("EM algorithm did not converge");
  } else {
    Rprintf("PL converged after %u iterations \n", it);
  }

  arma::vec pl = LL(beta, gamma, lambda, L, R, t, tind, X, S, 
                    knots, constantVE);

  return pl;
}
