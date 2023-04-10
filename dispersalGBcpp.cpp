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
int N_Eq_Str(CharacterVector ToCheck, CharacterVector Crit) { // check number cells vector matching character string 
  int n_equal = 0;
  for(int i = 0; i<ToCheck.size(); i++){
    if(ToCheck(i) == Crit(0)){ 
      n_equal++;
    }
  }
  return n_equal;
}

// [[Rcpp::export]]
IntegerVector sortInt(IntegerVector x) { // sort vector of integer
  IntegerVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}

// [[Rcpp::export]]
IntegerVector IntOrderIndex(IntegerVector x) {
  IntegerVector x_uniq = unique(x);
  IntegerVector x_uniqSorted = sortInt(x_uniq);
  IntegerVector Index(x.size());
  int p = 0; 
  for(int i = 0; i<x_uniqSorted.size(); i++){
    for(int j = 0; j<x.size(); j++){
      if(x(j) == x_uniqSorted(i)){
        Index(p) = j;
        p++;
      }
    }
  }
  return Index;
}


/////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List dispersalGB(// DataFrame, NumericVector
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
      IntegerVector CellsDisp_who(nDispLeft * 9);
      // hab freq indiv
      IntegerVector HabFreqDisp_barrier(nDispLeft);
      IntegerVector HabFreqDisp_matrix(nDispLeft);
      IntegerVector HabFreqDisp_disprep(nDispLeft);
      IntegerVector HabFreqDisp_who(nDispLeft);
      // vector recording whether next step is in matrix
      IntegerVector Mat_chosen(nDispLeft);
      for(int i = 0; i<nDispLeft; i++){
        for(int l2 = 0; l2<3; l2++){
          for(int l1 = 0; l1<3; l1++){
            CellsDisp_who(i * 9 + l1 + l2 * 3) = dispersers_who(dispersing(i));
            CellsDisp_x(i * 9 + l1 + l2 * 3) = dispersers_lastDispX(dispersing(i)) - 1 + l1;// because want square around indiv
            CellsDisp_y(i * 9 + l1 + l2 * 3) = dispersers_lastDispY(dispersing(i)) - 1 + l1;
            CellsDisp_ind(i * 9 + l1 + l2 * 3) = i;
            CellsDisp_hab(i * 9 + l1 + l2 * 3) = HabitatMap(CellsDisp_y(i * 9 + l1 + l2 * 3),CellsDisp_x(i * 9 + l1 + l2 * 3));
            // here count habitat occurences
            if(CellsDisp_hab(i * 9 + l1 + l2 * 3) == 0){ // barrier
              HabFreqDisp_who(i) = CellsDisp_who(i * 9 + l1 + l2 * 3);
              HabFreqDisp_barrier(i)++;
            }
            if(CellsDisp_hab(i * 9 + l1 + l2 * 3) == 2){ // matrix
              HabFreqDisp_who(i) = CellsDisp_who(i * 9 + l1 + l2 * 3);
              HabFreqDisp_matrix(i)++;
            }
            if(CellsDisp_hab(i * 9 + l1 + l2 * 3) == 3){ // dispersing
              HabFreqDisp_who(i) = CellsDisp_who(i * 9 + l1 + l2 * 3);
              HabFreqDisp_disprep(i)++;
            }
            if(CellsDisp_hab(i * 9 + l1 + l2 * 3) == 4){ // breeding
              HabFreqDisp_who(i) = CellsDisp_who(i * 9 + l1 + l2 * 3);
              HabFreqDisp_disprep(i)++;
            }
          }
        }
        //does indiv pick matrix type habitat for next move (pMat * freq hab 2)
        Mat_chosen(i) = R::rbinom(int_1, pMat * HabFreqDisp_matrix(i)); // number of samples is argument 2, ie matrix type
        //correction in case only habitat 2, otherwise indiv may be trown away
        if(HabFreqDisp_matrix(i) == 9){
          Mat_chosen(i) = 1;
        }
      }// end i loop indiv
      
      // find final number cell for dispersing individuals, taking matrix or disprep only 
      int nCellsDispLeft = 0;
      for(int nc = 0; nc<CellsDisp_x.size(); nc++){ // nc number of cells
        if((Mat_chosen(CellsDisp_ind(nc)) == 1) & (CellsDisp_hab(nc) == 2)){
          nCellsDispLeft++;
        }
        if((Mat_chosen(CellsDisp_ind(nc)) == 0) & ((CellsDisp_hab(nc) == 3) | (CellsDisp_hab(nc) == 4))){
          nCellsDispLeft++;
        }
      }
      
      // loop to pick only the right habit for further dispersion based on whether mat_chosen = 1
      IntegerVector nextCellsType_x(nCellsDispLeft);
      IntegerVector nextCellsType_y(nCellsDispLeft);
      IntegerVector nextCellsType_ind(nCellsDispLeft);
      IntegerVector nextCellsType_hab(nCellsDispLeft);
      IntegerVector nextCellsType_who(nCellsDispLeft);
      for(int i = 0, p = 0; i<CellsDisp_ind.size(); i++){
        if((Mat_chosen(CellsDisp_ind(i)) == 1) & (CellsDisp_hab(i) == 2)){
          //stop("Went into loop");
          nextCellsType_x(p) = CellsDisp_x(i);
          nextCellsType_y(p) = CellsDisp_y(i);
          nextCellsType_ind(p) = CellsDisp_ind(i);
          nextCellsType_hab(p) = CellsDisp_hab(i);
          nextCellsType_who(p) = CellsDisp_who(i);
          p++;
        }
        if((Mat_chosen(CellsDisp_ind(i)) == 0) & ((CellsDisp_hab(i) == 3) | (CellsDisp_hab(i) == 4))){
          nextCellsType_x(p) = CellsDisp_x(i);
          nextCellsType_y(p) = CellsDisp_y(i);
          nextCellsType_ind(p) = CellsDisp_ind(i);
          nextCellsType_hab(p) = CellsDisp_hab(i);
          nextCellsType_who(p) = CellsDisp_who(i);
          p++;
        }
      }// end second ind loop
      
      
      // part on potential correlation in movements, none on  first move but then some/////////////////////////////////////////
      if(step == 0){
        IntegerVector randLines_move = sample((nextCellsType_hab.size()), (nextCellsType_hab.size()), false) - 1; // number of samples is argument 2;
        IntegerVector ChosenCell_x(nDispLeft);
        IntegerVector ChosenCell_y(nDispLeft);
        IntegerVector ChosenCell_hab(nDispLeft);
        IntegerVector ChosenCell_ind(nDispLeft);
        IntegerVector ChosenCell_who(nDispLeft);
        for(int ind = 0; ind < nDispLeft; ind++){
          double p = 0.5;
          while(p<1){
            for(int l = 0; l<nextCellsType_hab.size(); l++){
              if(nextCellsType_ind(randLines_move(l)) == ind){ // here keeps last one, while was not working
                ChosenCell_x(ind) = nextCellsType_x(randLines_move(l));
                ChosenCell_y(ind) = nextCellsType_y(randLines_move(l));
                ChosenCell_hab(ind) = nextCellsType_hab(randLines_move(l));
                ChosenCell_ind(ind) = nextCellsType_ind(randLines_move(l));
                ChosenCell_who(ind) = nextCellsType_who(randLines_move(l));
                p = p+1;
              }
            }
          }
        }
        // List L_return = List::create(Named("nCellsDispLeft") = nCellsDispLeft,
        //                              _["nDispLeft"] = nDispLeft,
        //                              //_["randLines_move"] = randLines_move,
        //                              _["ChosenCell_x"] = ChosenCell_x,
        //                              _["ChosenCell_ind"] = ChosenCell_ind);
        // return  L_return;
      } //Works up to here
      else{ // i.e if step more than first
        // sort by ind and add current position next to potential ones in nextCellsType bectors
        IntegerVector nextCellsType_indSortIndex = IntOrderIndex(nextCellsType_ind);
        IntegerVector nextCellsType_indSorted = nextCellsType_ind.size();
        IntegerVector nextCellsType_xSorted = nextCellsType_ind.size();
        IntegerVector nextCellsType_ySorted = nextCellsType_ind.size();
        IntegerVector nextCellsType_habSorted = nextCellsType_ind.size();
        IntegerVector nextCellsType_xCurSorted = nextCellsType_ind.size();
        IntegerVector nextCellsType_yCurSorted = nextCellsType_ind.size();
        IntegerVector nextCellsType_whoSorted = nextCellsType_who.size();
        for(int l = 0; l< nextCellsType_ind.size(); l++){
          nextCellsType_indSorted(l) = nextCellsType_ind(nextCellsType_indSortIndex(l));
          nextCellsType_xSorted(l) = nextCellsType_x(nextCellsType_indSortIndex(l));
          nextCellsType_ySorted(l) = nextCellsType_y(nextCellsType_indSortIndex(l));
          nextCellsType_habSorted(l) = nextCellsType_hab(nextCellsType_indSortIndex(l));
          nextCellsType_whoSorted(l) = nextCellsType_who(nextCellsType_indSortIndex(l));
          for(int j = 0; j<dispersers_who.size(); j++){
            if(dispersers_who(j) == nextCellsType_whoSorted(l)){
              nextCellsType_xCurSorted(l) = dispersers_lastDispX(j);
              nextCellsType_yCurSorted(l) = dispersers_lastDispY(j);
            }
          }
        }
        List L_return = List::create(Named("nextCellsType_indSorted") =  nextCellsType_indSorted,
                                       _["nextCellsType_yCurSorted"] = nextCellsType_yCurSorted,
                                       _["nextCellsType_ySorted"] = nextCellsType_ySorted,
                                       _["nextCellsType_xCurSorted"] = nextCellsType_xCurSorted,
                                       _["nextCellsType_xSorted"] = nextCellsType_xSorted,
                                       _["nextCellsType_whoSorted"] = nextCellsType_whoSorted);
        return  L_return;
        
      }
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
