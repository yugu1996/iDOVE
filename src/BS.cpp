#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <ctime>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::rowvec BS(int t, arma::ivec knots, bool constantVE) {
  // define a row arma::vector {B_1(t), ..., B_K(t)}, based on the knots 
  // the slope of the last piece is 0 if constantVE = True
  int npc;
  double temp;
  arma::rowvec B;

  if (constantVE) {
    npc = knots.n_elem;
    temp = t>knots(npc-1) ? t-knots(npc-1) : 0;
  } else {
    npc = knots.n_elem+1;
    temp = 0;
  }

  B.set_size(npc);
  B(0) = t-temp;
  for (int i=1; i<npc; ++i) {
    B(i) = (t>knots(i-1) ? t-knots(i-1) : 0)-temp;
  }
  B = B*0.0329;
  return B;
}
