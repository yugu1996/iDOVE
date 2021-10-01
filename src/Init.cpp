#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <ctime>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

class Less {
  int _than;
public:
  Less(int th): _than(th){}
  bool operator()(int i){
    return i < _than;
  }
};

/* Returns {n x 3} matrix, containing the element of t
   nearest the entry, left interval, and right interval times

   L     in     {n} left interval time in days
   Rstar in     {n} right interval time in days modified to remove Inf
   W     in     {n} entry time in days
   t     in     {nTimes} time points for int
*/
// [[Rcpp::export]]
arma::imat Init(const arma::ivec& L, const arma::ivec& Rstar, 
                const arma::ivec & W, const arma::ivec& t) {

  arma::uvec Wsort = sort_index(W);
  arma::imat tind(L.n_elem,3); 

  for(auto i: Wsort) {
    tind(i,0) = std::count_if(t.begin(), t.end(), Less(W(i)));
    tind(i,1) = std::count_if(t.begin(), t.end(), Less(L(i)));
    tind(i,2) = std::count_if(t.begin(), t.end(), Less(Rstar(i)));
  }

  return tind;
}
