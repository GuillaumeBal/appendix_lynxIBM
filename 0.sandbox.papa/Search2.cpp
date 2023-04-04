#include <R.h>
//#include "Rcpp.h"
#include "RcppArmadillo.h"
using namespace arma; // Evite d'écrire arma::fonctionArmadillo
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]] // Importe la fonction qui suit dans l'environnement global de R


StringMatrix f_crea_matrice(
    StringVector Input_Vec1, 	//sexe
    LogicalVector Input_Vec2, 	//vivant
    StringVector Input_Vec3		//race
) {
  
  // Create an empty matrix with the desired dimensions
  StringMatrix disp(Input_Vec1.length(), 3);
  StringMatrix dispFem(Input_Vec1.length(), 3);
  
  // Replace the columns of the matrix by the input vectors
  disp( _ , 0) = Input_Vec1;
  disp( _ , 1) = Input_Vec2;
  disp( _ , 2) = Input_Vec3;
  
  if((disp.nrow())!= 0){
    int numL= 0;
    for (int i = 0; i<Input_Vec1.length(); i++) {      /* loop over all rows  */
      if (disp(i,0) == "F"){                   /* select rows  */
        dispFem(numL, _) = disp(i, _);			/* add row to dis */
        numL++;
      }
    }
    StringMatrix dispFem2(numL-1, 3);
    dispFem2 = dispFem(Range(0, numL-1) , Range(0, 2)); // Reference to sub matrix
    return dispFem2;
  }
  else{
  }
  
}
