#include <R.h>
//#include "Rcpp.h"
#include "RcppArmadillo.h"
using namespace arma; // Evite d'écrire arma::fonctionArmadillo
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]] // Importe la fonction qui suit dans l'environnement global de R


Rcpp::DataFrame towardsGB(
    DataFrame Disp
) {
  
  // some attributes of input dataframe
  //int Disp_list_size = Disp.size();
  
  // dataframe columns to vectors
  NumericVector Disp_xcor = Disp["xcor"], Disp_ycor = Disp["ycor"], Disp_who = Disp["who"], Disp_heading = Disp["heading"], Disp_prevX = Disp["prevX"], Disp_prevY = Disp["prevY"], Disp_age = Disp["age"], Disp_steps = Disp["steps"], Disp_lastDispX = Disp["lastDispX"], Disp_lastDispY = Disp["lastDispY"], Disp_nMat = Disp["nMat"], Disp_maleID = Disp["maleID"], Disp_nFem = Disp["nFem"], Disp_rdMortTerr = Disp["rdMortTerr"]; 
  CharacterVector Disp_breed = Disp["breed"], Disp_color = Disp["color"], Disp_pop = Disp["pop"], Disp_sex = Disp["sex"], Disp_status = Disp["status"];
  
  // compute number females in disp
  int n_fem_Disp = 0;
  for(int l = 0; l<Disp_sex.size(); l++) {      /* loop over all rows  */
    if(Disp_sex(l) == "F"){ 
      n_fem_Disp++;
    }
  }
  // if some females do find lines and create subset matrix
  if(n_fem_Disp > 1){
    IntegerVector index_fem_Disp = n_fem_Disp;
    for(int l = 0, j = 0; l<Disp_sex.size(); l++){
      if(Disp_sex(l) == "F"){
        index_fem_Disp[j] = l;
        j++;
      }
    }
    //NumericVector index_fem_Disp_rand_order = RcppArmadillo::sample(index_fem_Disp, index_fem_Disp.size(), FALSE, 1/index_fem_Disp.size()) ;
    IntegerVector index_fem_Disp_rand = sample(index_fem_Disp, index_fem_Disp.size());
    NumericVector DispFem_xcor = Disp_xcor[index_fem_Disp], DispFem_ycor = Disp_ycor[index_fem_Disp], DispFem_who = Disp_who[index_fem_Disp], DispFem_heading = Disp_heading[index_fem_Disp], DispFem_prevX = Disp_prevX[index_fem_Disp], DispFem_prevY = Disp_prevY[index_fem_Disp], DispFem_age = Disp_age[index_fem_Disp], DispFem_steps = Disp_steps[index_fem_Disp], DispFem_lastDispX = Disp_lastDispX[index_fem_Disp], DispFem_lastDispY = Disp_lastDispY[index_fem_Disp], DispFem_nMat = Disp_nMat[index_fem_Disp], DispFem_maleID = Disp_maleID[index_fem_Disp], DispFem_nFem = Disp_nFem[index_fem_Disp], DispFem_rdMortTerr = Disp_rdMortTerr[index_fem_Disp]; 
    CharacterVector DispFem_breed = Disp_breed[index_fem_Disp], DispFem_color = Disp_color[index_fem_Disp], DispFem_pop = Disp_pop[index_fem_Disp], DispFem_sex = Disp_sex[index_fem_Disp], DispFem_status = Disp_status[index_fem_Disp];
    DataFrame df_out = DataFrame::create(Named("index_fem_Disp") = index_fem_Disp , _["index_fem_Disp_rand"] = index_fem_Disp_rand);
    return df_out;
  }
  
  // 
  stop("'n_fem_disp' must be higher than 0.");
  
}
