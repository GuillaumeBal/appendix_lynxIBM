#include <R.h>
//#include "Rcpp.h"
#include "RcppArmadillo.h"
using namespace arma; // Evite d'écrire arma::fonctionArmadillo
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]] // Importe la fonction qui suit dans l'environnement global de R


Rcpp::NumericVector towardsGB( // DataFrame, NumericVector
    DataFrame Disp,
    IntegerMatrix HabitatMap,
    IntegerMatrix TerrMap
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
    IntegerVector index_fem_Disp_rand = sample(index_fem_Disp, index_fem_Disp.size());
    NumericVector DispFem_xcor = Disp_xcor[index_fem_Disp_rand], DispFem_ycor = Disp_ycor[index_fem_Disp_rand], DispFem_who = Disp_who[index_fem_Disp_rand], DispFem_heading = Disp_heading[index_fem_Disp_rand], DispFem_prevX = Disp_prevX[index_fem_Disp_rand], DispFem_prevY = Disp_prevY[index_fem_Disp_rand], DispFem_age = Disp_age[index_fem_Disp_rand], DispFem_steps = Disp_steps[index_fem_Disp_rand], DispFem_lastDispX = Disp_lastDispX[index_fem_Disp_rand], DispFem_lastDispY = Disp_lastDispY[index_fem_Disp_rand], DispFem_nMat = Disp_nMat[index_fem_Disp_rand], DispFem_maleID = Disp_maleID[index_fem_Disp_rand], DispFem_nFem = Disp_nFem[index_fem_Disp_rand], DispFem_rdMortTerr = Disp_rdMortTerr[index_fem_Disp_rand]; 
    CharacterVector DispFem_breed = Disp_breed[index_fem_Disp_rand], DispFem_color = Disp_color[index_fem_Disp_rand], DispFem_pop = Disp_pop[index_fem_Disp_rand], DispFem_sex = Disp_sex[index_fem_Disp_rand], DispFem_status = Disp_status[index_fem_Disp_rand];
    // a check of the sample fonction good functioning outouting the df
    //DataFrame df_out = DataFrame::create(Named("index_fem_Disp") = index_fem_Disp , _["index_fem_Disp_rand"] = index_fem_Disp_rand);
    DataFrame DispFem = DataFrame::create(Named("xcor") = DispFem_xcor, _["ycor"] = DispFem_ycor, _["who"] = DispFem_who, _["heading"] = DispFem_heading, _["prevX"] = DispFem_prevX, _["prevY"] = DispFem_prevY, _["breed"] = DispFem_breed, _["color"] = DispFem_color, _["pop"] = DispFem_pop, _["sex"] = DispFem_sex, _["age"] = DispFem_age, _["status"] = DispFem_status, _["steps"] = DispFem_steps, _["lastDispX"] = DispFem_lastDispX, _["lastDispY"] = DispFem_lastDispY, _["nMat"] = DispFem_nMat, _["maleID"] = DispFem_maleID, _["nFem"] = DispFem_nFem, _["rdMortTerr"] = DispFem_rdMortTerr);
    //return DispFem;
    
    // loop for female territorry search
    NumericVector vec_out = n_fem_Disp;
    for(int f = 0; f<n_fem_Disp; f++){
      if(HabitatMap(DispFem_lastDispX(f), DispFem_lastDispY(f)) == 4 & R_IsNA(TerrMap(DispFem_lastDispX(f), DispFem_lastDispX(f)))){
        //if(HabitatMap(1, 1) == 4 & R_IsNA(TerrMap(1, 1))){
        vec_out(f) = 0;
      }
      else{
        vec_out(f) = 1;
      }
    }
    return vec_out;
  }// end of if some dispersing females
  
  // 
  stop("'n_fem_disp' must be higher than 0.");
  
}
