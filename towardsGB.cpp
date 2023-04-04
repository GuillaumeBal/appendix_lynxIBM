#include <R.h>
//#include "Rcpp.h"
#include "RcppArmadillo.h"
using namespace arma; // Evite d'écrire arma::fonctionArmadillo
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]] // Importe la fonction qui suit dans l'environnement global de R


Rcpp::DataFrame towardsGB(
    DataFrame disp
) {
  
  // compute number females
  int n_fem = 0;
  CharacterVector disp_sex = disp["sex"];
  for(int i = 0; i<disp.nrow(); i++) {      /* loop over all rows  */
    if(disp_sex(i) == "F"){ 
      n_fem++;
    }
  }
  return n_fem
}
