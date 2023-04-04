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
  
  // some attributes of input dataframe
  int disp_list_size = disp.size();
  
  // compute number females in disp
  int n_fem_disp = 0;
  CharacterVector disp_sex = disp["sex"];
  for(int l = 0; l<disp_sex.size(); l++) {      /* loop over all rows  */
    if(disp_sex(l) == "F"){ 
      n_fem_disp++;
    }
  }
  // if some females do find lines and create subset matrix
  if(n_fem_disp > 1){
    IntegerVector index_disp_fem = n_fem_disp;
    for(int l = 0, j = 0; l<disp_sex.size(); l++){
      if(disp_sex(l) == "F"){
        index_disp_fem[j] = l;
        j++;
      }
    }
    List DispFem = disp[index_disp_fem];
    //for(int c = 0; c<disp_list_size; c++){
      //disp(1);
    //}
    return DispFem;
  }
  
  // 
  
  stop("'n_fem_disp' must be a positive value.");
}
