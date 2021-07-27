#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <ctime>
#include "BS.h"

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double Cox_noX_LL(const arma::vec& gamma, const arma::ivec& time, 
                  const arma::ivec& delta, const arma::ivec& W, 
                  const arma::ivec& S, const arma::ivec& knots, bool constantVE) {
  int nfail, j, n = delta.n_elem; 
  double tempt, S0, xbeta, ll = 0;
  arma::rowvec tmpe;

  for (int i=n-1; i>=0; ) {

    nfail = 0; 
    tempt = time(i);
    S0 = 0.0; 
    j = n-1; 

    while (j>=0 && time(j)>=tempt) {
      if (tempt>=W(j)) {
        if (S(j)<tempt) {
          tmpe = BS(tempt-S(j), knots, constantVE);
          xbeta = arma::as_scalar(tmpe*gamma);
        } else {
          xbeta = 0.0;
        }
        S0 += exp(xbeta);

        if (time(j)==tempt) {
          i--;
          if (delta(j)==1) {
            nfail++;
            ll += xbeta;
          }
        }
      }
      j--;
    }

    ll -= nfail*log(S0);

  }
  return ll;
}

void getUI_noX(arma::vec& U, arma::mat& I,
               const arma::vec& gamma, const arma::ivec& time, 
               const arma::ivec& delta, const arma::ivec& W, 
               const arma::ivec& S, const arma::ivec& knots, bool constantVE) {
  int nfail, j, n = delta.n_elem, p2 = gamma.n_elem; 
  double tempt, S0, xbeta, exbeta;
  arma::vec S1;
  arma::rowvec tmpe;
  arma::mat S2;

  for (int i=n-1; i>=0; ) {
    nfail = 0; 
    tempt = time(i);
    S0 = 0.0; 
    S1.zeros(p2); 
    S2.zeros(p2, p2);
    j = n-1; 
    while (j>=0 && time(j)>=tempt) {
      if (tempt>=W(j)) {
        if (S(j)<tempt) {
          tmpe = BS(tempt-S(j), knots, constantVE);
          xbeta = arma::as_scalar(tmpe*gamma);
          exbeta = exp(xbeta);
          S0 += exbeta;
          S1 += exbeta*tmpe.t();
          S2 += exbeta*tmpe.t()*tmpe;
        } else {
          tmpe.zeros(p2);
          exbeta = 1.0;
          S0 += exbeta;
        }

        if (time(j)==tempt) {
          i--;
          if (delta(j)==1) {
            nfail++;
            U += tmpe.t();
          }
        }
      }
      j--;
    }

    U -= nfail*S1/S0;
    I += nfail*(S2/S0-S1*S1.t()/pow(S0,2));
  }
}

arma::mat getw_noX(const arma::vec& gamma, 
                   const arma::ivec& time, const arma::ivec& delta, 
                   const arma::ivec& W, 
                   const arma::ivec& S, const arma::ivec& knots, bool constantVE) {

  // compute w matrix for robust variance estimator
  int j, m, n = delta.n_elem, p2 = gamma.n_elem; 
  double exbeta;
  arma::vec S0, D;
  arma::ivec t;
  arma::rowvec Zi;
  arma::mat S1, w;

  t = unique(time);
  m = t.n_elem; 
  D.zeros(m);
  S0.zeros(m); 
  S1.zeros(m, p2); 
  w.zeros(n, p2);
  for (int k=m-1; k>=0; --k) {
    j = n-1; 
    while (j>=0 && time(j)>=t(k)) {
      if (t(k)>=W(j)) {
        if (S(j)<t(k)) {
          Zi = BS(t(k)-S(j), knots, constantVE);
          exbeta = exp(arma::as_scalar(Zi*gamma));
          S1.row(k) += exbeta*Zi;
        } else {
          Zi.zeros(p2);
          exbeta = 1.0;
        }
        S0(k) += exbeta;

        if (time(j)==t(k)) {
          if (delta(j)==1) {
            D(k) += 1.0;
            w.row(j) += Zi;
          }
        }
      }
      j--;
    } // end j
  } // end k

  for (int i=0; i<n; ++i) {
    j = 0;
    while (j<m && time(i)>=t(j)) {
      if (t(j)>=W(i)) {
        if (S(i)<t(j)) {
          Zi = BS(t(j)-S(i), knots, constantVE);
          exbeta = exp(arma::as_scalar(Zi*gamma));
          w.row(i) -= Zi*exbeta/n/S0(j)*D(j);
        } else {
          exbeta = 1.0;
        }
        w.row(i) += S1.row(j)*exbeta/n/pow(S0(j),2)*D(j);

        if(t(j)==time(i) && delta(i)==1) {
          w.row(i) -= S1.row(j)/S0(j);
        }
      }
      j++;
    }
  }

  return w;
}

// [[Rcpp::export]]
List Cox_noX(arma::vec& gamma, const arma::ivec& time, const arma::ivec& delta, 
             const arma::ivec& W, const arma::ivec& S, 
             const arma::ivec& knots, 
             bool constantVE, double threshold=10^(-4), int maxit=5000) {
  // Standard Cox model with piecewise log-linear hazard ratio r(t) with 5 pieces
  int n = time.n_elem, p2 = gamma.n_elem, it = 0;
  double dif = 1.0, ll;
  arma::vec gamma0, U, temp;
  arma::rowvec XZ;
  arma::mat I, I_inv, w, covar(p2, p2);
  bool flag; // indicator whether I is invertible

  // begin: Newton-Raphson

  while (it<maxit && dif>threshold) {

    R_CheckUserInterrupt();

    it++;
    gamma0 = gamma;
    U.zeros(p2); 
    I.zeros(p2, p2);

    getUI_noX(U, I, gamma, time, delta, W, S, knots, constantVE);

    flag = pinv(I_inv, I);

    if (!flag) {
      stop("Newton-Raphson terminated due to singular information matrix");
    }

    temp = I_inv*U;
    gamma += temp;
    dif = norm(gamma-gamma0, 2);
  }
  // end: Newton-Raphson

  // compute robust variance estimator

  ll = Cox_noX_LL(gamma, time, delta, W, S, knots, constantVE);
  w = getw_noX(gamma, time, delta, W, S, knots, constantVE);
  covar = I_inv*w.t()*w*I_inv;

  List to_return(2); 
  // to_return: ll, covariance matrix
  to_return[0] = ll;
  to_return[1] = covar;
  return to_return;
}
