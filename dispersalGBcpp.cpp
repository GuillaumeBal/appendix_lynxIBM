#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

////////////////////////////////////////////////////////////////////////////////
// series of helper function for lynx script

// [[Rcpp::export]]
int N_Eq_Str(CharacterVector ToCheck, CharacterVector Crit) {
  int n_equal = 0;
  for(int i = 0; i<ToCheck.size(); i++){
    if(ToCheck(i) == Crit(0)){ 
      n_equal++;
    }
  }
  return n_equal;
}



// [[Rcpp::export]]
IntegerVector dispersalGB(// DataFrame, NumericVector
    DataFrame lynx,
    int sMaxPs, // dispersal
    IntegerMatrix HabitatMap,
    double pMat
){
  
  ///////////////////////////////////////////////////////////////////////////////////////
  // TRICKS AND FUNCTIONS FOR CODE
  
  // some integer def because IntegerVector v(1) = 1 does not work
  int int_0 = 0;
  int int_1 = 1;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  // LYNX DATA PROCESSING
  
  // attach all lynx data
  NumericVector lynx_xcor = lynx["xcor"], lynx_ycor = lynx["ycor"], lynx_who = lynx["who"], lynx_heading = lynx["heading"], lynx_prevX = lynx["prevX"], lynx_prevY = lynx["prevY"], lynx_age = lynx["age"], lynx_lastDispX = lynx["lastDispX"], lynx_lastDispY = lynx["lastDispY"], lynx_nMat = lynx["nMat"], lynx_maleID = lynx["maleID"], lynx_nFem = lynx["nFem"], lynx_rdMortTerr = lynx["rdMortTerr"]; 
  CharacterVector lynx_breed = lynx["breed"], lynx_color = lynx["color"], lynx_pop = lynx["pop"], lynx_sex = lynx["sex"], lynx_status = lynx["status"];
  
  // create disperser and non disperser var /////////////////////////////////////////////////////////////////////
  
  //define some dim
  int nLynx = lynx_xcor.size();
  int nDisp = N_Eq_Str(lynx_status, "disp");
  // define dispersers vectors
  IntegerVector dispersers_xcor(nDisp), dispersers_ycor(nDisp), dispersers_who(nDisp), dispersers_heading(nDisp), dispersers_prevX(nDisp), dispersers_prevY(nDisp), dispersers_age(nDisp), dispersers_lastDispX(nDisp), dispersers_lastDispY(nDisp), dispersers_nMat(nDisp), dispersers_maleID(nDisp), dispersers_nFem(nDisp), dispersers_rdMortTerr(nDisp); 
  CharacterVector dispersers_breed(nDisp), dispersers_color(nDisp), dispersers_pop(nDisp), dispersers_sex(nDisp), dispersers_status(nDisp);
  // define residents vectors
  IntegerVector residents_xcor(nLynx - nDisp), residents_ycor(nLynx - nDisp), residents_who(nLynx - nDisp), residents_heading(nLynx - nDisp), residents_prevX(nLynx - nDisp), residents_prevY(nLynx - nDisp), residents_age(nLynx - nDisp), residents_lastDispX(nLynx - nDisp), residents_lastDispY(nLynx - nDisp), residents_nMat(nLynx - nDisp), residents_maleID(nLynx - nDisp), residents_nFem(nLynx - nDisp), residents_rdMortTerr(nLynx - nDisp); 
  CharacterVector residents_breed(nLynx - nDisp), residents_color(nLynx - nDisp), residents_pop(nLynx - nDisp), residents_sex(nLynx - nDisp), residents_status(nLynx - nDisp);
  // make loop to fill residetns and dispersers vectors
  for(int i = 0, i_disp = 0, i_ndisp = 0 ; i<lynx_status.size(); i++){
    if(lynx_status(i) == "disp"){
      dispersers_xcor(i_disp) = lynx_xcor(i), dispersers_ycor(i_disp) = lynx_ycor(i), dispersers_who(i_disp) = lynx_who(i), dispersers_heading(i_disp) = lynx_heading(i), dispersers_prevX(i_disp) = lynx_prevX(i), dispersers_prevY(i_disp) = lynx_prevY(i), dispersers_breed(i_disp) = lynx_breed(i), dispersers_color(i_disp) = lynx_color(i), dispersers_pop(i_disp) = lynx_pop(i), dispersers_sex(i_disp) = lynx_sex(i), dispersers_age(i_disp) = lynx_age(i), dispersers_status(i_disp) = lynx_status(i), dispersers_lastDispX(i_disp) = lynx_lastDispX(i), dispersers_lastDispY(i_disp) = lynx_lastDispY(i), dispersers_nMat(i_disp) = lynx_nMat(i), dispersers_maleID(i_disp) = lynx_maleID(i), dispersers_nFem(i_disp) = lynx_nFem(i), dispersers_rdMortTerr(i_disp) = lynx_rdMortTerr(i);
      i_disp++;
    }
    else{
      residents_xcor(i_disp) = lynx_xcor(i), residents_ycor(i_disp) = lynx_ycor(i), residents_who(i_disp) = lynx_who(i), residents_heading(i_disp) = lynx_heading(i), residents_prevX(i_disp) = lynx_prevX(i), residents_prevY(i_disp) = lynx_prevY(i), residents_breed(i_disp) = lynx_breed(i), residents_color(i_disp) = lynx_color(i), residents_pop(i_disp) = lynx_pop(i), residents_sex(i_disp) = lynx_sex(i), residents_age(i_disp) = lynx_age(i), residents_status(i_disp) = lynx_status(i), residents_lastDispX(i_disp) = lynx_lastDispX(i), residents_lastDispY(i_disp) = lynx_lastDispY(i), residents_nMat(i_disp) = lynx_nMat(i), residents_maleID(i_disp) = lynx_maleID(i), residents_nFem(i_disp) = lynx_nFem(i), residents_rdMortTerr(i_disp) = lynx_rdMortTerr(i);
      i_ndisp ++;
    }
  }
  //return dispersers_xcor;
  
  // work on disperser from now on //////////////////////////////////////////////////////////
  
  IntegerVector stepsDisp = sample(sMaxPs, nDisp, true); // number of samples is argument 2
  int maxDisp = max(stepsDisp);
  // find indiv that are still dispersing, ie have some more steps to do withing the day
  for(int step = 0; step<maxDisp; step++){
    IntegerVector dispersing(nDisp);
    if(nDisp>0){
      for(int d = 0, dp = 0; d<dispersers_status.size(); d++){
        if(step<stepsDisp(d)){
          dispersing(dp) = d;
          dp++;
        }
      }
    }
    //return dispersing; //works up to here without crashing after a while
    int nDispLeft = dispersing.size(); // i.e are there some steps left for the indiv ?
    if(nDispLeft>0){
      //mat with X / Y / id ind / Habitat cell
      IntegerVector CellsDisp_x(nDispLeft * 9);
      IntegerVector CellsDisp_y(nDispLeft * 9);
      IntegerVector CellsDisp_ind(nDispLeft * 9);
      IntegerVector CellsDisp_hab(nDispLeft * 9);
      // hab freq indiv
      IntegerVector HabFreqDisp_barrier(nDispLeft);
      IntegerVector HabFreqDisp_matrix(nDispLeft);
      IntegerVector HabFreqDisp_disprep(nDispLeft);
      // vector recording whether next step is in matrix
      IntegerVector Mat_chosen(nDispLeft);
      for(int i = 0; i<nDispLeft; i++){
        for(int l2 = 0; l2<3; l2++){
          for(int l1 = 0; l1<3; l1++){
            CellsDisp_x(i * 9 + l1 + l2 * 3) = dispersers_lastDispX(dispersing(i)) - 1 + l1;// because want square around indiv
            CellsDisp_y(i * 9 + l1 + l2 * 3) = dispersers_lastDispY(dispersing(i)) - 1 + l1;
            CellsDisp_ind(i * 9 + l1 + l2 * 3) = i;
            CellsDisp_hab(i * 9 + l1 + l2 * 3) = HabitatMap(CellsDisp_x(i * 9 + l1 + l2 * 3), CellsDisp_y(i * 9 + l1 + l2 * 3));
            // here count habitat occurences
            if(CellsDisp_hab(i * 9 + l1 + l2 * 3) == 0){ // barrier
              HabFreqDisp_barrier(i)++;
            }
            if(CellsDisp_hab(i * 9 + l1 + l2 * 3) == 2){ // matrix
              HabFreqDisp_matrix(i)++;
            }
            if(CellsDisp_hab(i * 9 + l1 + l2 * 3) == 3){ // dispersing
              HabFreqDisp_disprep(i)++;
            }
            if(CellsDisp_hab(i * 9 + l1 + l2 * 3) == 4){ // breeding
              HabFreqDisp_disprep(i)++;
            }
          }
        }
        //does indiv pick matrix type habitat for next move (pMat * freq hab 2)
        Mat_chosen(i) = R::rbinom(int_1, pMat * HabFreqDisp_matrix(i)); // number of samples is argument 2, ie matrix type
      }// end i loop indiv
      //return CellsDisp_hab; // works up to here
      
      
      // find final number cell for dispersing individuals
      int nCellsDispLefts = 0;
      for(int nc = 0; nc<CellsDisp_x.size(); nc++){ // nc number of cells
        if((Mat_chosen(CellsDisp_ind(nc)) == 1) & (CellsDisp_hab(nc) == 2)){
          nCellsDispLefts++;
        }
        if((Mat_chosen(CellsDisp_ind(nc)) == 0) & ((CellsDisp_hab(nc) == 3) | (CellsDisp_hab(nc) == 4))){
          nCellsDispLefts++;
        }
      }
      return nCellsDispLefts; 
      
      // loop to pick only the right habit for further dispersion based on whether mat_chosen = 1
      for(int i = 0, p = 0; i<nDispLeft; i++){
        // cell left to disperse
        IntegerVector nextCellsType_x(nCellsDispLefts);
        IntegerVector nextCellsType_y(nCellsDispLefts);
        IntegerVector nextCellsType_ind(nCellsDispLefts);
        IntegerVector nextCellsType_hab(nCellsDispLefts);
        for(int nci = 0; nci<9; nci++){// nci number cell per indiv, fill final dispersal matrix base on choice
          if((Mat_chosen[i] == 1) & (CellsDisp_hab(nci + i * nci) == 2)){
            nextCellsType_x(p) = CellsDisp_x(nci + i * nci);
            nextCellsType_y(p) = CellsDisp_y(nci + i * nci);
            nextCellsType_ind(p) = CellsDisp_ind(nci + i * nci);
            nextCellsType_hab(p) = CellsDisp_hab(nci + i * nci);
            p++;
          }else{
            if((Mat_chosen[i] == 0) & (CellsDisp_hab(nci + i * nci) == 3)){
              //stop("got within");
              nextCellsType_x(p) = CellsDisp_x(nci + i * nci);
              nextCellsType_y(p) = CellsDisp_y(nci + i * nci);
              nextCellsType_ind(p) = CellsDisp_ind(nci + i * nci);
              nextCellsType_hab(p) = CellsDisp_hab(nci + i * nci);
              p++;
            }else{
              if((Mat_chosen[i] == 0) & (CellsDisp_hab(nci + i * nci) == 4)){
                //stop("got within");
                nextCellsType_x(p) = CellsDisp_x(nci + i * nci);
                nextCellsType_y(p) = CellsDisp_y(nci + i * nci);
                nextCellsType_ind(p) = CellsDisp_ind(nci + i * nci);
                nextCellsType_hab(p) = CellsDisp_hab(nci + i * nci);
                p++;
              }
            }
          }
          return nextCellsType_hab;
        }
      }// end second ind loop
    }// of if dispersers left
  }// end step loop
} // end of function

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

//*** R
//timesTwo(42)
//*/
