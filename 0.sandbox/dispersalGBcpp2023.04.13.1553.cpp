#include <Rcpp.h>
using namespace Rcpp;
#include <math.h>       /* atan2 */
#define PI 3.14159265
#include <algorithm>

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

// [[Rcpp::export]]
IntegerVector WhichAbove(IntegerVector ToCheck, int Crit) { // check number cells vector matching character string 
  int n_above = 0;
  for(int i = 0; i<ToCheck.size();i++){
    if(ToCheck(i)>Crit){
      n_above++;
    }
  }
  IntegerVector which_vec(n_above);
  int p = 0;
  for(int i = 0; i<ToCheck.size(); i++){
    if(ToCheck(i)>Crit){
      which_vec(p) = i;
      p++;
    }
  }
  return which_vec;
}

// [[Rcpp::export]]
IntegerVector WhichEqual(IntegerVector ToCheck, int Crit) { // // give index cells valua above crit for IntegerVector  
  int n_equal = 0;
  for(int i = 0; i<ToCheck.size();i++){
    if(ToCheck(i)==Crit){
      n_equal++;
    }
  }
  IntegerVector which_vec(n_equal);
  int p = 0;
  for(int i = 0; i<ToCheck.size(); i++){
    if(ToCheck(i)==Crit){
      which_vec(p) = i;
      p++;
    }
  }
  return which_vec;
}

// [[Rcpp::export]]
IntegerVector IntVecSubIndex(IntegerVector ToSub, IntegerVector PosToKeep) { // subset integer vector based on index
  IntegerVector kept(PosToKeep.size());
  for(int i = 0; i<PosToKeep.size(); i++){
    kept(i) = ToSub(PosToKeep(i));
  }
  return kept;
}

// [[Rcpp::export]]
IntegerVector WhichInSetInt(IntegerVector ToCheck, IntegerVector subset) { // check cells within subset for integer 
  int n_in = 0;
  for(int i = 0; i<ToCheck.size(); i++){
    for(int j = 0; j<subset.size(); j++){
      if(ToCheck(i)==subset(j)){
        n_in++;
      }
    }
  }
  IntegerVector which_vec(n_in);
  int p = 0;
  for(int i = 0; i<ToCheck.size(); i++){
    for(int j = 0; j<subset.size(); j++){
      if(ToCheck(i)==subset(j)){
        which_vec(p) = i;
        p++;
      }
    }
  }
  return which_vec;
}

// [[Rcpp::export]]
// here randomly shuffling line before picking one per unique number of x
IntegerVector IntPosOneOfEach(IntegerVector x){
  IntegerVector randLines_move = sample(x.size(), x.size(), false) - 1; // number of samples is argument 2; -1 to make it start at 0
  IntegerVector unique_x = unique(x);
  IntegerVector unique_x_ordered = clone(unique_x);
  unique_x_ordered.sort(false);
  IntegerVector Chosen_Line(unique_x.size());
  for(int ind = 0; ind < unique_x_ordered.size(); ind++){
    double p = 0.5;
    while(p<1){
      for(int l = 0; l<x.size(); l++){
        if(x(randLines_move(l)) == unique_x_ordered(ind)){ // here keeps last one, while was not working
          Chosen_Line(ind) = randLines_move(l);
          p = p + 1;
        }
      }
    }
  }
  return Chosen_Line;
}

// [[Rcpp::export]]
// towards simple,  convert radians to degrees
int towards_simple_unique(int x_cur, int y_cur, int x_to, int y_to){
  int degrees = 0;
  degrees = atan2(x_to - x_cur, y_to - y_cur) * (180 / PI);
  if(degrees<0){
    degrees = degrees + 360;
  }
  return degrees;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List dispersalGB(// DataFrame, NumericVector
    DataFrame lynx,
    int sMaxPs, // dispersal
    IntegerMatrix HabitatMap,
    double pMat,
    double pCorr,
    int nMatMax,
    IntegerMatrix connectivityMap,
    NumericMatrix roadMortMap,
    int corrFactorDisp,
    int floorTimeSim,
    int startSimYear,
    IntegerVector ncoll_ncoll,
    IntegerVector ncoll_time,
    DataFrame deadLynxColl
){
  
  ///////////////////////////////////////////////////////////////////////////////////////
  // TRICKS AND FUNCTIONS FOR CODE
  
  // some integer def because IntegerVector v(1) = 1 does not work
  int int_0 = 0;
  int int_1 = 1;
  int int_2 = 2;
  int int_3 = 3;
  int int_4 = 4;
  int int_5 = 5;
  int int_100 = 100;
  
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
  
  IntegerVector dispersers_steps = sample(sMaxPs, nDisp, true); // number of samples is argument 2
  int maxDisp = max(dispersers_steps);
  // find indiv that are still dispersing, ie have some more steps to do withing the day
  for(int step = 0; step<maxDisp; step++){
    IntegerVector dispersing(nDisp);
    if(nDisp>0){
      for(int d = 0, dp = 0; d<dispersers_status.size(); d++){
        if(step<dispersers_steps(d)){
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
      IntegerVector CellsDisp_steps(nDispLeft * 9);
      IntegerVector CellsDisp_nMat(nDispLeft * 9);
      IntegerVector CellsDisp_heading(nDispLeft * 9);
      // hab freq indiv
      IntegerVector HabFreqDisp_barrier(nDispLeft);
      IntegerVector HabFreqDisp_matrix(nDispLeft);
      IntegerVector HabFreqDisp_disprep(nDispLeft);
      IntegerVector HabFreqDisp_who(nDispLeft);
      IntegerVector HabFreqDisp_steps(nDispLeft);
      IntegerVector HabFreqDisp_nMat(nDispLeft);
      IntegerVector HabFreqDisp_heading(nDispLeft);
      // vector recording whether next step is in matrix
      IntegerVector Mat_Chosen(nDispLeft);
      for(int i = 0; i<nDispLeft; i++){
        for(int l2 = 0; l2<3; l2++){
          for(int l1 = 0; l1<3; l1++){
            CellsDisp_who(i * 9 + l1 + l2 * 3) = dispersers_who(dispersing(i));
            CellsDisp_steps(i * 9 + l1 + l2 * 3) = dispersers_steps(dispersing(i));
            CellsDisp_nMat(i * 9 + l1 + l2 * 3) = dispersers_nMat(dispersing(i));
            CellsDisp_heading(i * 9 + l1 + l2 * 3) = dispersers_heading(dispersing(i));
            CellsDisp_x(i * 9 + l1 + l2 * 3) = dispersers_lastDispX(dispersing(i)) - 1 + l1;// because want square around indiv
            CellsDisp_y(i * 9 + l1 + l2 * 3) = dispersers_lastDispY(dispersing(i)) - 1 + l1;
            CellsDisp_ind(i * 9 + l1 + l2 * 3) = i;
            CellsDisp_hab(i * 9 + l1 + l2 * 3) = HabitatMap(CellsDisp_y(i * 9 + l1 + l2 * 3),CellsDisp_x(i * 9 + l1 + l2 * 3));
            // here count habitat occurences
            if(CellsDisp_hab(i * 9 + l1 + l2 * 3) == 0){ // barrier
              HabFreqDisp_who(i) = CellsDisp_who(i * 9 + l1 + l2 * 3);
              HabFreqDisp_steps(i) = CellsDisp_steps(i * 9 + l1 + l2 * 3);
              HabFreqDisp_nMat(i) = CellsDisp_nMat(i * 9 + l1 + l2 * 3);
              HabFreqDisp_heading(i) = CellsDisp_nMat(i * 9 + l1 + l2 * 3);
              HabFreqDisp_barrier(i)++;
            }
            if(CellsDisp_hab(i * 9 + l1 + l2 * 3) == 2){ // matrix
              HabFreqDisp_who(i) = CellsDisp_who(i * 9 + l1 + l2 * 3);
              HabFreqDisp_steps(i) = CellsDisp_steps(i * 9 + l1 + l2 * 3);
              HabFreqDisp_nMat(i) = CellsDisp_nMat(i * 9 + l1 + l2 * 3);
              HabFreqDisp_heading(i) = CellsDisp_heading(i * 9 + l1 + l2 * 3);
              HabFreqDisp_matrix(i)++;
            }
            if(CellsDisp_hab(i * 9 + l1 + l2 * 3) == 3){ // dispersing
              HabFreqDisp_who(i) = CellsDisp_who(i * 9 + l1 + l2 * 3);
              HabFreqDisp_steps(i) = CellsDisp_steps(i * 9 + l1 + l2 * 3);
              HabFreqDisp_nMat(i) = CellsDisp_nMat(i * 9 + l1 + l2 * 3);
              HabFreqDisp_heading(i) = CellsDisp_nMat(i * 9 + l1 + l2 * 3);
              HabFreqDisp_disprep(i)++;
            }
            if(CellsDisp_hab(i * 9 + l1 + l2 * 3) == 4){ // breeding
              HabFreqDisp_who(i) = CellsDisp_who(i * 9 + l1 + l2 * 3);
              HabFreqDisp_steps(i) = CellsDisp_steps(i * 9 + l1 + l2 * 3);
              HabFreqDisp_nMat(i) = CellsDisp_nMat(i * 9 + l1 + l2 * 3);
              HabFreqDisp_heading(i) = CellsDisp_nMat(i * 9 + l1 + l2 * 3);
              HabFreqDisp_disprep(i)++;
            }
          }
        }
        //does indiv pick matrix type habitat for next move (pMat * freq hab 2)
        Mat_Chosen(i) = R::rbinom(int_1, pMat * HabFreqDisp_matrix(i)); // number of samples is argument 2, ie matrix type
        //correction in case only habitat 2, otherwise indiv may be trown away
        if(HabFreqDisp_matrix(i) == 9){
          Mat_Chosen(i) = 1;
        }
      }// end i loop indiv
      
      // find final number cell for dispersing individuals, taking matrix or disprep only 
      int nCellsDispLeft = 0;
      for(int nc = 0; nc<CellsDisp_x.size(); nc++){ // nc number of cells
        if((Mat_Chosen(CellsDisp_ind(nc)) == 1) & (CellsDisp_hab(nc) == 2)){
          nCellsDispLeft++;
        }
        if((Mat_Chosen(CellsDisp_ind(nc)) == 0) & ((CellsDisp_hab(nc) == 3) | (CellsDisp_hab(nc) == 4))){
          nCellsDispLeft++;
        }
      }
      
      // loop to pick only the right habit for further dispersion based on whether mat_Chosen = 1
      IntegerVector nextCellsType_x(nCellsDispLeft);
      IntegerVector nextCellsType_y(nCellsDispLeft);
      IntegerVector nextCellsType_ind(nCellsDispLeft);
      IntegerVector nextCellsType_hab(nCellsDispLeft);
      IntegerVector nextCellsType_who(nCellsDispLeft);
      IntegerVector nextCellsType_steps(nCellsDispLeft);
      IntegerVector nextCellsType_nMat(nCellsDispLeft);
      IntegerVector nextCellsType_heading(nCellsDispLeft);
      for(int i = 0, p = 0; i<CellsDisp_ind.size(); i++){
        if((Mat_Chosen(CellsDisp_ind(i)) == 1) & (CellsDisp_hab(i) == 2)){
          //stop("Went into loop");
          nextCellsType_x(p) = CellsDisp_x(i);
          nextCellsType_y(p) = CellsDisp_y(i);
          nextCellsType_ind(p) = CellsDisp_ind(i);
          nextCellsType_hab(p) = CellsDisp_hab(i);
          nextCellsType_who(p) = CellsDisp_who(i);
          nextCellsType_steps(p) = CellsDisp_steps(i);
          nextCellsType_nMat(p) = CellsDisp_nMat(i);
          nextCellsType_heading(p) = CellsDisp_heading(i);
          p++;
        }
        if((Mat_Chosen(CellsDisp_ind(i)) == 0) & ((CellsDisp_hab(i) == 3) | (CellsDisp_hab(i) == 4))){
          nextCellsType_x(p) = CellsDisp_x(i);
          nextCellsType_y(p) = CellsDisp_y(i);
          nextCellsType_ind(p) = CellsDisp_ind(i);
          nextCellsType_hab(p) = CellsDisp_hab(i);
          nextCellsType_who(p) = CellsDisp_who(i);
          nextCellsType_steps(p) = CellsDisp_steps(i);
          nextCellsType_nMat(p) = CellsDisp_nMat(i);
          nextCellsType_heading(p) = CellsDisp_heading(i);
          p++;
        }
      }// end second ind loop
      
      // final steps elements def here fopr push back
      IntegerVector ChosenCells_x(0);
      IntegerVector ChosenCells_y(0);
      IntegerVector ChosenCells_xCur(0);
      IntegerVector ChosenCells_yCur(0);
      IntegerVector ChosenCells_hab(0);
      IntegerVector ChosenCells_ind(0);
      IntegerVector ChosenCells_who(0);
      IntegerVector ChosenCells_steps(0);
      IntegerVector ChosenCells_IsMoveCorr(0);
      IntegerVector ChosenCells_nMat(0);
      IntegerVector ChosenCells_heading(0);
      
      // part on potential correlation in movements, none on first move but then some/////////////////////////////////////////
      
      if(step == 0){
        // here randomly shuffling line before picking one per ind
        IntegerVector randLines_move = sample((nextCellsType_hab.size()), (nextCellsType_hab.size()), false) - 1; // number of samples is argument 2;
        for(int ind = 0; ind < nDispLeft; ind++){
          double p = 0.5;
          while(p<1){
            for(int l = 0; l<nextCellsType_hab.size(); l++){
              if(nextCellsType_ind(randLines_move(l)) == ind){ // here keeps last one, while was not working
                ChosenCells_x.push_back(nextCellsType_x(randLines_move(l)));
                ChosenCells_y.push_back(nextCellsType_y(randLines_move(l)));
                ChosenCells_hab.push_back(nextCellsType_hab(randLines_move(l)));
                ChosenCells_ind.push_back(nextCellsType_ind(randLines_move(l)));
                ChosenCells_who.push_back(nextCellsType_who(randLines_move(l)));
                ChosenCells_steps.push_back(nextCellsType_steps(randLines_move(l)));
                ChosenCells_nMat.push_back(nextCellsType_nMat(randLines_move(l)));
                ChosenCells_heading.push_back(nextCellsType_heading(randLines_move(l)));
                ChosenCells_IsMoveCorr.push_back(0);
                p = p + 1;
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
      
      else{ // i.e if step more than first //////////////////////////////////////////////////
        
        //Have to define here and push_back values to be able to use latter on outside loop where it is filled up
        // final move of indiv with no correlated movements
        IntegerVector ChosenCellsNoCorr_ind(0);
        IntegerVector ChosenCellsNoCorr_x(0);
        IntegerVector ChosenCellsNoCorr_y(0);
        IntegerVector ChosenCellsNoCorr_hab(0);
        IntegerVector ChosenCellsNoCorr_xCur(0);
        IntegerVector ChosenCellsNoCorr_yCur(0);
        IntegerVector ChosenCellsNoCorr_who(0);
        IntegerVector ChosenCellsNoCorr_steps(0);
        IntegerVector ChosenCellsNoCorr_nMat(0);
        IntegerVector ChosenCellsNoCorr_heading(0);
        IntegerVector ChosenCellsNoCorr_IsMoveCorr(0);
        // final moce infiv with correlatedmovement
        IntegerVector ChosenCellsYesCorr_who(0);
        IntegerVector ChosenCellsYesCorr_prefDir(0);
        IntegerVector ChosenCellsYesCorr_ind(0);
        IntegerVector ChosenCellsYesCorr_x(0);
        IntegerVector ChosenCellsYesCorr_y(0);
        IntegerVector ChosenCellsYesCorr_hab(0);
        IntegerVector ChosenCellsYesCorr_xCur(0);
        IntegerVector ChosenCellsYesCorr_yCur(0);
        IntegerVector ChosenCellsYesCorr_steps(0);
        IntegerVector ChosenCellsYesCorr_nMat(0);
        IntegerVector ChosenCellsYesCorr_heading(0);
        IntegerVector ChosenCellsYesCorr_IsMoveCorr(0);
        
        // sort by ind and add current position next to potential ones in nextCellsType vectors
        IntegerVector nextCellsType_indSortIndex = IntOrderIndex(nextCellsType_ind);
        IntegerVector nextCellsType_indSorted = nextCellsType_ind.size();
        IntegerVector nextCellsType_xSorted = nextCellsType_ind.size();
        IntegerVector nextCellsType_ySorted = nextCellsType_ind.size();
        IntegerVector nextCellsType_habSorted = nextCellsType_ind.size();
        IntegerVector nextCellsType_xCurSorted = nextCellsType_ind.size();
        IntegerVector nextCellsType_yCurSorted = nextCellsType_ind.size();
        IntegerVector nextCellsType_whoSorted = nextCellsType_who.size();
        IntegerVector nextCellsType_stepsSorted = nextCellsType_steps.size();
        IntegerVector nextCellsType_nMatSorted = nextCellsType_nMat.size();
        IntegerVector nextCellsType_headingSorted = nextCellsType_heading.size();
        for(int l = 0; l< nextCellsType_ind.size(); l++){
          nextCellsType_indSorted(l) = nextCellsType_ind(nextCellsType_indSortIndex(l));
          nextCellsType_xSorted(l) = nextCellsType_x(nextCellsType_indSortIndex(l));
          nextCellsType_ySorted(l) = nextCellsType_y(nextCellsType_indSortIndex(l));
          nextCellsType_habSorted(l) = nextCellsType_hab(nextCellsType_indSortIndex(l));
          nextCellsType_whoSorted(l) = nextCellsType_who(nextCellsType_indSortIndex(l));
          nextCellsType_stepsSorted(l) = nextCellsType_steps(nextCellsType_indSortIndex(l));
          nextCellsType_nMatSorted(l) = nextCellsType_nMat(nextCellsType_indSortIndex(l));
          nextCellsType_headingSorted(l) = nextCellsType_heading(nextCellsType_indSortIndex(l));
          for(int j = 0; j<dispersers_who.size(); j++){
            if(dispersers_who(j) == nextCellsType_whoSorted(l)){
              nextCellsType_xCurSorted(l) = dispersers_lastDispX(j);
              nextCellsType_yCurSorted(l) = dispersers_lastDispY(j);
            }
          }
        }
        // dispersers with some steps left, ie final matrix
        IntegerVector WStepsLeft = WhichAbove(nextCellsType_stepsSorted, step - 1);//- 1 because I want >= behavior from > function 
        IntegerVector nextCellsType_indF = IntVecSubIndex(nextCellsType_indSorted, WStepsLeft);
        IntegerVector nextCellsType_xF = IntVecSubIndex( nextCellsType_xSorted, WStepsLeft);
        IntegerVector nextCellsType_yF = IntVecSubIndex(nextCellsType_ySorted, WStepsLeft);
        IntegerVector nextCellsType_habF = IntVecSubIndex(nextCellsType_habSorted, WStepsLeft);
        IntegerVector nextCellsType_xCurF = IntVecSubIndex(nextCellsType_xCurSorted, WStepsLeft);
        IntegerVector nextCellsType_yCurF = IntVecSubIndex(nextCellsType_yCurSorted, WStepsLeft);
        IntegerVector nextCellsType_whoF = IntVecSubIndex(nextCellsType_whoSorted, WStepsLeft);
        IntegerVector nextCellsType_stepsF = IntVecSubIndex(nextCellsType_stepsSorted, WStepsLeft);
        IntegerVector nextCellsType_nMatF = IntVecSubIndex(nextCellsType_nMatSorted, WStepsLeft);
        IntegerVector nextCellsType_headingF = IntVecSubIndex(nextCellsType_headingSorted, WStepsLeft);
        
        //individual with correlated movement ?
        IntegerVector IsMoveCorr_who = unique(nextCellsType_whoF);
        int nDispF = IsMoveCorr_who.size();
        IntegerVector IsMoveCorr(nDispF);
        int nCorr1 = 0;
        for(int i = 0; i<nDispF; i++){
          IsMoveCorr(i) = R::rbinom(int_1, pCorr); // not a prob but a boolean 0/1
          if(IsMoveCorr(i) == int_1){
            nCorr1++;
          }
        }
        // assign to table for easier work
        IntegerVector nextCellsType_IsMoveCorrF(nextCellsType_stepsF.size());
        for(int i= 0; i<nDispF; i++){
          for(int j = 0; j<nextCellsType_stepsF.size(); j++){
            if(nextCellsType_whoF(j) == IsMoveCorr_who(i)){
              nextCellsType_IsMoveCorrF(j) = IsMoveCorr(i);
            }
          }
        }
        
        if(nCorr1< nDispF){ // i.e if some indiv with uncorrelated movement
          //stop("got in loop");
          IntegerVector noCorr_Lind = WhichEqual(nextCellsType_IsMoveCorrF, int_0); 
          int nCorr0 = noCorr_Lind.size();
          // define vector to fill
          IntegerVector nextCellsTypeNoCorr_indF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_xF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_yF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_habF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_xCurF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_yCurF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_whoF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_stepsF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_nMatF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_headingF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_IsMoveCorrF(nCorr0);
          for(int l = 0; l<nCorr0;l++){
            nextCellsTypeNoCorr_indF(l) = nextCellsType_indF(noCorr_Lind(l));
            nextCellsTypeNoCorr_xF(l) = nextCellsType_xF(noCorr_Lind(l));
            nextCellsTypeNoCorr_yF(l) = nextCellsType_yF(noCorr_Lind(l));
            nextCellsTypeNoCorr_habF(l) = nextCellsType_habF(noCorr_Lind(l));
            nextCellsTypeNoCorr_xCurF(l) = nextCellsType_xCurF(noCorr_Lind(l));
            nextCellsTypeNoCorr_yCurF(l) = nextCellsType_yCurF(noCorr_Lind(l));
            nextCellsTypeNoCorr_whoF(l) = nextCellsType_whoF(noCorr_Lind(l));
            nextCellsTypeNoCorr_stepsF(l) = nextCellsType_stepsF(noCorr_Lind(l));
            nextCellsTypeNoCorr_nMatF(l) = nextCellsType_nMatF(noCorr_Lind(l));
            nextCellsTypeNoCorr_headingF(l) = nextCellsType_headingF(noCorr_Lind(l));
            nextCellsTypeNoCorr_IsMoveCorrF(l) = nextCellsType_IsMoveCorrF(noCorr_Lind(l));
          }
          // now pick just one cell to move to for each ind
          IntegerVector UniqueLinesIndNoCorr = IntPosOneOfEach(nextCellsTypeNoCorr_indF);
          //IntegerVector ChosenCellsNoCorr_ind(UniqueLinesIndNoCorr.size());
          // IntegerVector ChosenCellsNoCorr_x(UniqueLinesIndNoCorr.size());
          // IntegerVector ChosenCellsNoCorr_y(UniqueLinesIndNoCorr.size());
          // IntegerVector ChosenCellsNoCorr_hab(UniqueLinesIndNoCorr.size());
          // IntegerVector ChosenCellsNoCorr_xCur(UniqueLinesIndNoCorr.size());
          // IntegerVector ChosenCellsNoCorr_yCur(UniqueLinesIndNoCorr.size());
          // IntegerVector ChosenCellsNoCorr_who(UniqueLinesIndNoCorr.size());
          // IntegerVector ChosenCellsNoCorr_steps(UniqueLinesIndNoCorr.size());
          // IntegerVector ChosenCellsNoCorr_nMat(UniqueLinesIndNoCorr.size());
          // IntegerVector ChosenCellsNoCorr_heading(UniqueLinesIndNoCorr.size());
          // IntegerVector ChosenCellsNoCorr_IsMoveCorr(UniqueLinesIndNoCorr.size());
          for(int i = 0; i<UniqueLinesIndNoCorr.size(); i++){
            ChosenCellsNoCorr_ind.push_back(nextCellsTypeNoCorr_indF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_x.push_back(nextCellsTypeNoCorr_xF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_y.push_back(nextCellsTypeNoCorr_yF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_hab.push_back(nextCellsTypeNoCorr_habF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_xCur.push_back(nextCellsTypeNoCorr_xCurF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_yCur.push_back(nextCellsTypeNoCorr_yCurF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_who.push_back(nextCellsTypeNoCorr_whoF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_steps.push_back(nextCellsTypeNoCorr_stepsF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_nMat.push_back(nextCellsTypeNoCorr_nMatF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_heading.push_back(nextCellsTypeNoCorr_headingF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_IsMoveCorr.push_back(nextCellsTypeNoCorr_IsMoveCorrF(UniqueLinesIndNoCorr(i)));
          }
          // List L_return = List::create(Named("ChosenCellsNoCorr_who") = ChosenCellsNoCorr_who,
          //                              _["ChosenCellsNoCorr_ind"] = ChosenCellsNoCorr_ind,
          //                              _["ChosenCellsNoCorr_yCur"] = ChosenCellsNoCorr_yCur);
          //return  L_return;
        }// end noCorr move indiv, works up to here
        
        // now deal with indiv with correlated movement
        if(nCorr1 > int_0){
          //stop("got into loop");
          IntegerVector YesCorr_Lind = WhichEqual(nextCellsType_IsMoveCorrF, int_1);
          //define vector to fill
          IntegerVector nextCellsTypeYesCorr_indF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_xF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_yF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_habF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_xCurF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_yCurF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_whoF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_stepsF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_nMatF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_headingF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_IsMoveCorrF(YesCorr_Lind.size());
          for(int l = 0; l<YesCorr_Lind.size(); l++){
            nextCellsTypeYesCorr_indF(l) = nextCellsType_indF(YesCorr_Lind(l));
            nextCellsTypeYesCorr_xF(l) = nextCellsType_xF(YesCorr_Lind(l));
            nextCellsTypeYesCorr_yF(l) = nextCellsType_yF(YesCorr_Lind(l));
            nextCellsTypeYesCorr_habF(l) = nextCellsType_habF(YesCorr_Lind(l));
            nextCellsTypeYesCorr_xCurF(l) = nextCellsType_xCurF(YesCorr_Lind(l));
            nextCellsTypeYesCorr_yCurF(l) = nextCellsType_yCurF(YesCorr_Lind(l));
            nextCellsTypeYesCorr_whoF(l) = nextCellsType_whoF(YesCorr_Lind(l));
            nextCellsTypeYesCorr_stepsF(l) = nextCellsType_stepsF(YesCorr_Lind(l));
            nextCellsTypeYesCorr_nMatF(l) = nextCellsType_nMatF(YesCorr_Lind(l));
            nextCellsTypeYesCorr_headingF(l) = nextCellsType_headingF(YesCorr_Lind(l));
            nextCellsTypeYesCorr_IsMoveCorrF(l) = nextCellsType_IsMoveCorrF(YesCorr_Lind(l));
          }
          // add dir to nextCellsTypeYesCorr
          IntegerVector nextCellsTypeYesCorr_DirF(nextCellsTypeYesCorr_IsMoveCorrF.size());
          IntegerVector nextCellsTypeYesCorr_prefDirF(nextCellsTypeYesCorr_IsMoveCorrF.size());
          for(int l = 0; l<nextCellsTypeYesCorr_IsMoveCorrF.size(); l++){
            nextCellsTypeYesCorr_DirF(l) = towards_simple_unique(nextCellsTypeYesCorr_xCurF(l), nextCellsTypeYesCorr_yCurF(l), 
                                      nextCellsTypeYesCorr_xF(l), nextCellsTypeYesCorr_yF(l));
            nextCellsTypeYesCorr_prefDirF(l) = int_1; // set to one then modify
            if((nextCellsTypeYesCorr_DirF(l) >= 45) & (nextCellsTypeYesCorr_DirF(l) <= 315)){
              nextCellsTypeYesCorr_prefDirF(l) = int_2;
            }
            if((nextCellsTypeYesCorr_DirF(l) >= 90) & (nextCellsTypeYesCorr_DirF(l) <= 270)){
              nextCellsTypeYesCorr_prefDirF(l) = int_3;
            }
            if((nextCellsTypeYesCorr_DirF(l) >= 135) & (nextCellsTypeYesCorr_DirF(l) <= 225)){
              nextCellsTypeYesCorr_prefDirF(l) = int_4;
            }
            if(nextCellsTypeYesCorr_DirF(l) == 180){
              nextCellsTypeYesCorr_prefDirF(l) = int_5;
            }
            if((nextCellsTypeYesCorr_xCurF(l) == nextCellsTypeYesCorr_xF(l)) & (nextCellsTypeYesCorr_yCurF(l) == nextCellsTypeYesCorr_yF(l))){
              nextCellsTypeYesCorr_prefDirF(l) = int_3;
            }
          }
          // keep only one line per indiv, with lowest value of rank
          IntegerVector unique_nextCellsTypeYesCorr_whoF = unique(nextCellsTypeYesCorr_whoF);
          unique_nextCellsTypeYesCorr_whoF = sortInt(unique_nextCellsTypeYesCorr_whoF);
          // IntegerVector ChosenCellsYesCorr_who(unique_nextCellsTypeYesCorr_whoF.size());
          // IntegerVector ChosenCellsYesCorr_prefDir(unique_nextCellsTypeYesCorr_whoF.size());
          // IntegerVector ChosenCellsYesCorr_ind(unique_nextCellsTypeYesCorr_whoF.size());
          // IntegerVector ChosenCellsYesCorr_x(unique_nextCellsTypeYesCorr_whoF.size());
          // IntegerVector ChosenCellsYesCorr_y(unique_nextCellsTypeYesCorr_whoF.size());
          // IntegerVector ChosenCellsYesCorr_hab(unique_nextCellsTypeYesCorr_whoF.size());
          // IntegerVector ChosenCellsYesCorr_xCur(unique_nextCellsTypeYesCorr_whoF.size());
          // IntegerVector ChosenCellsYesCorr_yCur(unique_nextCellsTypeYesCorr_whoF.size());
          // IntegerVector ChosenCellsYesCorr_steps(unique_nextCellsTypeYesCorr_whoF.size());
          // IntegerVector ChosenCellsYesCorr_nMat(unique_nextCellsTypeYesCorr_nMatF.size());
          // IntegerVector ChosenCellsYesCorr_heading(unique_nextCellsTypeYesCorr_headingF.size());
          // IntegerVector ChosenCellsYesCorr_IsMoveCorr(unique_nextCellsTypeYesCorr_whoF.size());
          // subset for one of the lower prefdir values
          for(int i = 0; i<unique_nextCellsTypeYesCorr_whoF.size(); i++){
            ChosenCellsYesCorr_prefDir.push_back(int_100);// set to 100 to be able to replace by values within table
            ChosenCellsYesCorr_who.push_back(int_100);
            ChosenCellsYesCorr_ind.push_back(int_100);
            ChosenCellsYesCorr_x.push_back(int_100);
            ChosenCellsYesCorr_y.push_back(int_100);
            ChosenCellsYesCorr_hab.push_back(int_100);
            ChosenCellsYesCorr_xCur.push_back(int_100);
            ChosenCellsYesCorr_yCur.push_back(int_100);
            ChosenCellsYesCorr_steps.push_back(int_100);
            ChosenCellsYesCorr_nMat.push_back(int_100);
            ChosenCellsYesCorr_heading.push_back(int_100);
            ChosenCellsYesCorr_IsMoveCorr.push_back(int_100);
            for(int l = 0; l<nextCellsTypeYesCorr_whoF.size() ; l++){
              if((nextCellsTypeYesCorr_prefDirF(l) < ChosenCellsYesCorr_prefDir(i)) &
                 (nextCellsTypeYesCorr_whoF(l) == unique_nextCellsTypeYesCorr_whoF(i))){
                ChosenCellsYesCorr_prefDir(i) = nextCellsTypeYesCorr_prefDirF(l);
                ChosenCellsYesCorr_who(i) = nextCellsTypeYesCorr_whoF(l);
                ChosenCellsYesCorr_ind(i) = nextCellsTypeYesCorr_indF(l);
                ChosenCellsYesCorr_x(i) = nextCellsTypeYesCorr_xF(l);
                ChosenCellsYesCorr_y(i) = nextCellsTypeYesCorr_yF(l);
                ChosenCellsYesCorr_hab(i) = nextCellsTypeYesCorr_habF(l);
                ChosenCellsYesCorr_xCur(i) = nextCellsTypeYesCorr_xCurF(l);
                ChosenCellsYesCorr_yCur(i) = nextCellsTypeYesCorr_yCurF(l);
                ChosenCellsYesCorr_steps(i) = nextCellsTypeYesCorr_stepsF(l);
                ChosenCellsYesCorr_nMat(i) = nextCellsTypeYesCorr_nMatF(l);
                ChosenCellsYesCorr_heading(i) = nextCellsTypeYesCorr_headingF(l);
                ChosenCellsYesCorr_IsMoveCorr(i) = nextCellsTypeYesCorr_IsMoveCorrF(l);
              }
            }
          }
          
          // List L_return = List::create(Named("nextCellsType_IsMoveCorrF") = nextCellsType_IsMoveCorrF,
          //                              //_["YesCorr_Lind"] = YesCorr_Lind,
          //                              _["nCorr1"] = nCorr1,
          //                              _["nextCellsTypeYesCorr_indF"] = nextCellsTypeYesCorr_indF,
          //                              _["nextCellsTypeYesCorr_whoF"] = nextCellsTypeYesCorr_whoF,
          //                              _["nextCellsTypeYesCorr_DirF"] = nextCellsTypeYesCorr_DirF,
          //                              _["nextCellsTypeYesCorr_prefDirF"] = nextCellsTypeYesCorr_prefDirF,
          //                              _["ChosenCellsYesCorr_prefDir"] = ChosenCellsYesCorr_prefDir,
          //                              _["ChosenCellsYesCorr_who"] = ChosenCellsYesCorr_who,
          //                              _["ChosenCellsYesCorr_hab"] = ChosenCellsYesCorr_hab);
          // return  L_return;
          
          // get together moves for correlated and uncorrelated individuals
          int nLChosenCellsYesCorr = ChosenCellsYesCorr_who.size();
          int nLChosenCellsNoCorr = ChosenCellsNoCorr_who.size(); //ChosenCellsNoCorr
          int nChosenCells = nLChosenCellsYesCorr + nLChosenCellsNoCorr;
          // IntegerVector ChosenCells_who(nChosenCells);
          // IntegerVector ChosenCells_ind(nChosenCells);
          // IntegerVector ChosenCells_hab(nChosenCells);
          // IntegerVector ChosenCells_x(nChosenCells);
          // IntegerVector ChosenCells_y(nChosenCells);
          // IntegerVector ChosenCells_xCur(nChosenCells);
          // IntegerVector ChosenCells_yCur(nChosenCells);
          // IntegerVector ChosenCells_nMat(nChosenCells);
          // IntegerVector ChosenCells_heading(nChosenCells);
          for(int l = 0; l<nChosenCells; l++){
            if(l < nLChosenCellsYesCorr){
              ChosenCells_who.push_back(ChosenCellsYesCorr_who(l));
              ChosenCells_ind.push_back(ChosenCellsYesCorr_ind(l));
              ChosenCells_hab.push_back(ChosenCellsYesCorr_hab(l));
              ChosenCells_x.push_back(ChosenCellsYesCorr_x(l));
              ChosenCells_y.push_back(ChosenCellsYesCorr_y(l));
              ChosenCells_xCur.push_back(ChosenCellsYesCorr_xCur(l));
              ChosenCells_yCur.push_back(ChosenCellsYesCorr_yCur(l));
              ChosenCells_nMat.push_back(ChosenCellsYesCorr_nMat(l));
              ChosenCells_heading.push_back(ChosenCellsYesCorr_heading(l));
              ChosenCells_IsMoveCorr.push_back(ChosenCellsYesCorr_IsMoveCorr(l));
            }else{
              ChosenCells_who.push_back(ChosenCellsNoCorr_who(l - nLChosenCellsYesCorr));
              ChosenCells_ind.push_back(ChosenCellsNoCorr_ind(l - nLChosenCellsYesCorr));
              ChosenCells_hab.push_back(ChosenCellsNoCorr_hab(l - nLChosenCellsYesCorr));
              ChosenCells_x.push_back(ChosenCellsNoCorr_x(l - nLChosenCellsYesCorr));
              ChosenCells_y.push_back(ChosenCellsNoCorr_y(l - nLChosenCellsYesCorr));
              ChosenCells_xCur.push_back(ChosenCellsNoCorr_xCur(l - nLChosenCellsYesCorr));
              ChosenCells_yCur.push_back(ChosenCellsNoCorr_yCur(l - nLChosenCellsYesCorr));
              ChosenCells_nMat.push_back(ChosenCellsNoCorr_nMat(l - nLChosenCellsYesCorr));
              ChosenCells_heading.push_back(ChosenCellsNoCorr_heading(l - nLChosenCellsYesCorr));
              ChosenCells_IsMoveCorr.push_back(ChosenCellsNoCorr_IsMoveCorr(l - nLChosenCellsYesCorr));
            }
          }
        }
        
        // work on chosenMat and chosenDisp matrices and their processing /////////////////////////////////////////
        IntegerVector MatInd = WhichEqual(ChosenCells_hab, int_2);
        int nMatInd = MatInd.size();
        int nDispInd = ChosenCells_nMat.size() - nMatInd;
        IntegerVector ChosenMat_who(nMatInd);
        IntegerVector ChosenMat_ind(nMatInd);
        IntegerVector ChosenMat_hab(nMatInd);
        IntegerVector ChosenMat_x(nMatInd);
        IntegerVector ChosenMat_y(nMatInd);
        IntegerVector ChosenMat_xCur(nMatInd);
        IntegerVector ChosenMat_yCur(nMatInd);
        IntegerVector ChosenMat_nMat(nMatInd);
        IntegerVector ChosenMat_heading(nMatInd);
        // disp
        IntegerVector ChosenDisp_who(nDispInd);
        IntegerVector ChosenDisp_ind(nDispInd);
        IntegerVector ChosenDisp_hab(nDispInd);
        IntegerVector ChosenDisp_x(nDispInd);
        IntegerVector ChosenDisp_y(nDispInd);
        IntegerVector ChosenDisp_xCur(nDispInd);
        IntegerVector ChosenDisp_yCur(nDispInd);
        IntegerVector ChosenDisp_nMat(nDispInd);
        IntegerVector ChosenDisp_heading(nDispInd);
        for(int l = 0, p = 0, q = 0; l<ChosenCells_nMat.size(); l++){
          if(ChosenCells_hab(l) == 2){
            if(nMatInd>0){
              ChosenMat_who(p) = ChosenCells_who(l);
              ChosenMat_ind(p) = ChosenCells_ind(l);
              ChosenMat_hab(p) = ChosenCells_hab(l);
              ChosenMat_x(p) = ChosenCells_x(l);
              ChosenMat_y(p) = ChosenCells_y(l);
              ChosenMat_xCur(p) = ChosenCells_xCur(l);
              ChosenMat_yCur(p) = ChosenCells_yCur(l);
              ChosenMat_nMat(p) = ChosenCells_nMat(l) + 1;
              // code memory in movement bit here
              if((ChosenMat_nMat(p) + 1) == nMatMax){
                ChosenMat_x(p) = ChosenMat_xCur(p); // there is really that
                ChosenMat_y(p) = ChosenMat_yCur(p);
              }if(ChosenMat_nMat(p) == nMatMax){ // reset nMat
                ChosenMat_nMat(p) = 0;
                ChosenMat_heading(p) = ChosenMat_heading(p) + 180;
              }
              p++;
            }
          }
          if((ChosenCells_hab(l) == 3) | (ChosenCells_hab(l) == 4)){
            if(nDispInd>0){
              ChosenDisp_who(q) = ChosenCells_who(l);
              ChosenDisp_ind(q) = ChosenCells_ind(l);
              ChosenDisp_hab(q) = ChosenCells_hab(l);
              ChosenDisp_x(q) = ChosenCells_x(l);
              ChosenDisp_y(q) = ChosenCells_y(l);
              ChosenDisp_xCur(q) = ChosenCells_xCur(l);
              ChosenDisp_yCur(q) = ChosenCells_yCur(l);
              ChosenDisp_nMat(q) = 0;
              q++;
            }
          }
        }
        
        // process ChosenMat lynxMemory already integrated within loops above
        // same for pxcor /lastdisp update  
        
        // now reupdate chosen cells matrix by binding ChosenDisp and ChosenMat 
        for(int l = 0; l<ChosenCells_who.size(); l++){
          if(l<nDispInd){
            ChosenCells_who(l) = ChosenDisp_who(l);
            ChosenCells_ind(l) = ChosenDisp_ind(l);
            ChosenCells_hab(l) =ChosenDisp_hab(l);
            ChosenCells_x(l) = ChosenDisp_x(l);
            ChosenCells_y(l) = ChosenDisp_y(l);
            ChosenCells_xCur(l) = ChosenDisp_xCur(l);
            ChosenCells_yCur(l) = ChosenDisp_yCur(l);
            ChosenCells_nMat(l) = ChosenDisp_nMat(l);
          }
          else{
            ChosenCells_who(l) = ChosenMat_who(l - nDispInd);
            ChosenCells_ind(l) = ChosenMat_ind(l - nDispInd);
            ChosenCells_hab(l) =ChosenMat_hab(l - nDispInd);
            ChosenCells_x(l) = ChosenMat_x(l - nDispInd);
            ChosenCells_y(l) = ChosenMat_y(l - nDispInd);
            ChosenCells_xCur(l) = ChosenMat_xCur(l - nDispInd);
            ChosenCells_yCur(l) = ChosenMat_yCur(l - nDispInd);
            ChosenCells_nMat(l) = ChosenMat_nMat(l - nDispInd);
          }
        }
        // update connectivity map with +1 when dipserser stp on cell
        for(int l = 0; l<ChosenCells_x.size(); l++){
          connectivityMap(ChosenCells_y(l) , ChosenCells_x(l)) +=1; 
        }
        
        //////////////////////////////////////////////////////////////////////////////////////////////
        // bit on mortality
        IntegerVector deathRoad(ChosenCells_x.size());
        for(int l = 0; l<ChosenCells_x.size(); l++){
          deathRoad(l) = R::rbinom(int_1, (roadMortMap(ChosenCells_y(l), ChosenCells_x(l)) / corrFactorDisp));
          if(floorTimeSim == startSimYear){ // cannot die first year
            deathRoad(l) = int_0;
          } 
        }
        ncoll_ncoll.push_back(sum(deathRoad));
        ncoll_time.push_back(floorTimeSim);
        for(int l = 0; l<ChosenCells_x.size(); l++){
          if(roadMortMap(ChosenCells_y(l), ChosenCells_x(l)) == 1){ // force death on border to simulate emigration but does not count as death in line above
            deathRoad(l) = int_1;
          }
        }
        // save some data on dead individual
        for(int l = 0; l < deathRoad.size(); l++){
          if(deathRoad(l) == int_1){
            deadLynxColl.push_back(ChosenCells_who(l), "who");
          }
        }
        
        
        
        List L_return = List::create(Named("nextCellsType_indF") = nextCellsType_indF,
                                     _["nextCellsType_IsMoveCorrF"] = nextCellsType_IsMoveCorrF,
                                     _["nextCellsType_yF"] = nextCellsType_yF,
                                     _["IsMoveCorr"] = IsMoveCorr,
                                     _["nCorr1"] = nCorr1,
                                     _["ChosenCellsYesCorr_prefDir"] = ChosenCellsYesCorr_prefDir,
                                     _["ChosenCellsYesCorr_ind"] = ChosenCellsYesCorr_ind,
                                     _["ChosenCells_y"] = ChosenCells_y,
                                     _["ChosenCells_who"] = ChosenCells_who,
                                     _["ChosenCells_ind"] = ChosenCells_ind,
                                     _["ChosenCells_hab"] = ChosenCells_hab,
                                     _["ChosenCells_IsMoveCorr"] = ChosenCells_IsMoveCorr,
                                     _["MatInd"] = MatInd,
                                     _["deathRoad"] = deathRoad);
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
