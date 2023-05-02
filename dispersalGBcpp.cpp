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
// check number cells vector matching character string
int N_Eq_Str(StringVector ToCheck, StringVector Crit) {  
  int n_equal = 0;
  for(int i = 0; i<ToCheck.size(); i++){
    if(ToCheck(i) == Crit(0)){ 
      n_equal++;
    }
  }
  return n_equal;
}

// [[Rcpp::export]]
// sort vector of integer
IntegerVector sortInt(IntegerVector x) { 
  IntegerVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}

// [[Rcpp::export]]
// get vector of order to sort by increasing value
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
// found which values are above crit
IntegerVector WhichAbove(IntegerVector ToCheck, int Crit) {  
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
// found which value is equal to crit
IntegerVector WhichEqual(IntegerVector ToCheck, int Crit) { 
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
// subset integer vector based on index
IntegerVector IntVecSubIndex(IntegerVector ToSub, IntegerVector PosToKeep) { 
  IntegerVector kept(PosToKeep.size());
  for(int i = 0; i<PosToKeep.size(); i++){
    kept(i) = ToSub(PosToKeep(i));
  }
  return kept;
}

// [[Rcpp::export]]
// check cells within subset for integer 
IntegerVector WhichInSetInt(IntegerVector ToCheck, IntegerVector subset) { 
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

// [[Rcpp::export]]
// adjacent cells coordinates
List UniqFreeAdjCellsRandOrd(IntegerVector x_coords, IntegerVector y_coords, IntegerMatrix Matrix){
  IntegerVector deltaX = {-1, 0, 1};
  IntegerVector deltaY = {-1, 0, 1};
  int nColMat = Matrix.ncol();
  int nRowMat = Matrix.nrow();
  IntegerVector AdjX(0);
  IntegerVector AdjY(0);
  IntegerVector CellInd(0);
  int new_y;
  int new_x;
  int new_index;
  for(int z = 0; z<x_coords.size(); z++){
    for(int l = 0; l<deltaY.size(); l++){
      for(int c = 0; c<deltaX.size(); c++){
        new_y = y_coords(z) + deltaY(l);
        new_x = x_coords(z) + deltaX(c);
        new_index = (new_x + 1) * nRowMat + new_y;
        if((new_y>=0) & (new_y<nRowMat) & (new_x>=0) & (new_x<nColMat) &
           (Matrix(new_y, new_x) == 0) &
           ((deltaY(l) != 0) | (deltaX(c) != 0))){ // add to keep only unoccupied cell  s
          AdjX.push_back(new_x);
          AdjY.push_back(new_y);
          CellInd.push_back(new_index);
          //stop("Went into loop");
        }
      }
    }
  }
  // find unique index and then pick just one of each
  IntegerVector UniqCell = unique(CellInd);
  int nUniqCell = UniqCell.size();
  IntegerVector keptCells(nUniqCell);
  for(int i = 0; i<nUniqCell; i++){
    for(int j = 0; j<CellInd.size(); j++){
      if(CellInd(j) == UniqCell(i)){
        keptCells(i) = j;
        //break;
      }
    }
  }
  IntegerVector randorder = sample(AdjX.size(), AdjX.size(), false) - 1;
  IntegerVector keptCellsRandOrder = keptCells[randorder];
  IntegerVector AdjX_left = AdjX[keptCellsRandOrder];
  IntegerVector AdjY_left = AdjY[keptCellsRandOrder];
  IntegerVector CellInd_left = CellInd[keptCellsRandOrder];
  List L_return = List::create(Named("AdjX") = AdjX_left,
                               _["AdjY"] = AdjY_left,
                               _["CellInd"] = CellInd_left);
  // return L_return;
  return L_return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List dispersalGB(// DataFrame, NumericVector
    DataFrame lynx_r,
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
    DataFrame deadLynxColl,
    DataFrame deadDisp,
    IntegerMatrix TerrMap,
    IntegerMatrix availCellsUpdatedRas,
    IntegerMatrix popDist,
    int coreTerrSizeFAlps,
    int coreTerrSizeFJura,
    int coreTerrSizeFVosgesPalatinate,
    int coreTerrSizeFBlackForest,
    bool returnDistances,
    bool allowOverlap
){
  
  ///////////////////////////////////////////////////////////////////////////////////////
  // TRICKS AND FUNCTIONS FOR CODE
  
  DataFrame lynx = clone(lynx_r); // if no cloning, rccp keeps a link with table within R and resizing is problematic sometimes
  
  // some integer def because IntegerVector v(1) = 1 does not work
  int int_0 = 0;
  int int_1 = 1;
  int int_2 = 2;
  int int_3 = 3;
  int int_4 = 4;
  int int_5 = 5;
  int int_100 = 100;
  
  // defined to survey matrix filling dimension issue
  IntegerVector deathRoadOld(0);
  int nDispLeftOld(0);
  
  ///////////////////////////////////////////////////////////////////////////////////////
  // LYNX DATA PROCESSING
  
  int nLynx = lynx.nrow();
  int nLynInit = nLynx;
  // attach all lynx data
  IntegerVector lynx_xcor = lynx["xcor"], lynx_ycor = lynx["ycor"], lynx_who = lynx["who"], lynx_heading = lynx["heading"], lynx_prevX = lynx["prevX"], lynx_prevY = lynx["prevY"], lynx_age = lynx["age"], lynx_lastDispX = lynx["lastDispX"], lynx_lastDispY = lynx["lastDispY"], lynx_nMat = lynx["nMat"], lynx_maleID = lynx["maleID"], lynx_nFem = lynx["nFem"], lynx_rdMortTerr = lynx["rdMortTerr"]; 
  CharacterVector lynx_breed = lynx["breed"], lynx_color = lynx["color"], lynx_pop = lynx["pop"], lynx_sex = lynx["sex"], lynx_status = lynx["status"];
  // set a number of steps for all even though used only for dispersers
  IntegerVector lynx_steps = sample(sMaxPs - 1, nLynx, true) + 1; // number of samples is argument 2
  StringVector lynx_dispNow(nLynx);
  for(int l = 0; l<nLynx; l++){
    if((lynx_steps(l) > 0) & (lynx_status(l) == "disp")){
      lynx_dispNow(l) = "yes";
    }else{
      lynx_dispNow(l) = "no";
    }
  }
  
  // create disperser and non disperser var /////////////////////////////////////////////////////////////////////
  
  //define some dim
  int nDisp = N_Eq_Str(lynx_status, "disp");
  int nDispLeft = N_Eq_Str(lynx_dispNow, "yes");
  int nRes = N_Eq_Str(lynx_dispNow, "no");
  // define dispersers vectors
  IntegerVector dispersers_xcor(nDisp), dispersers_steps(nDisp), dispersers_ycor(nDisp), dispersers_who(nDisp), dispersers_heading(nDisp), dispersers_prevX(nDisp), dispersers_prevY(nDisp), dispersers_age(nDisp), dispersers_lastDispX(nDisp), dispersers_lastDispY(nDisp), dispersers_nMat(nDisp), dispersers_maleID(nDisp), dispersers_nFem(nDisp), dispersers_rdMortTerr(nDisp); 
  StringVector dispersers_breed(nDisp), dispersers_color(nDisp), dispersers_pop(nDisp), dispersers_sex(nDisp), dispersers_status(nDisp);
  // define residents vectors
  IntegerVector residents_xcor(nLynx - nDisp), residents_steps(nLynx - nDisp), residents_ycor(nLynx - nDisp), residents_who(nLynx - nDisp), residents_heading(nLynx - nDisp), residents_prevX(nLynx - nDisp), residents_prevY(nLynx - nDisp), residents_age(nLynx - nDisp), residents_lastDispX(nLynx - nDisp), residents_lastDispY(nLynx - nDisp), residents_nMat(nLynx - nDisp), residents_maleID(nLynx - nDisp), residents_nFem(nLynx - nDisp), residents_rdMortTerr(nLynx - nDisp); 
  StringVector residents_breed(nLynx - nDisp), residents_color(nLynx - nDisp), residents_pop(nLynx - nDisp), residents_sex(nLynx - nDisp), residents_status(nLynx - nDisp);
  // make loop to fill residetns and dispersers vectors
  for(int i = 0, i_disp = 0, i_ndisp = 0 ; i<lynx_status.size(); i++){
    if(lynx_dispNow(i) == "yes"){
      dispersers_xcor(i_disp) = lynx_xcor(i), dispersers_steps(i_disp) = lynx_steps(i), dispersers_ycor(i_disp) = lynx_ycor(i), dispersers_who(i_disp) = lynx_who(i), dispersers_heading(i_disp) = lynx_heading(i), dispersers_prevX(i_disp) = lynx_prevX(i), dispersers_prevY(i_disp) = lynx_prevY(i), dispersers_breed(i_disp) = lynx_breed(i), dispersers_color(i_disp) = lynx_color(i), dispersers_pop(i_disp) = lynx_pop(i), dispersers_sex(i_disp) = lynx_sex(i), dispersers_age(i_disp) = lynx_age(i), dispersers_status(i_disp) = lynx_status(i), dispersers_lastDispX(i_disp) = lynx_lastDispX(i), dispersers_lastDispY(i_disp) = lynx_lastDispY(i), dispersers_nMat(i_disp) = lynx_nMat(i), dispersers_maleID(i_disp) = lynx_maleID(i), dispersers_nFem(i_disp) = lynx_nFem(i), dispersers_rdMortTerr(i_disp) = lynx_rdMortTerr(i);
      i_disp++;
    }
    else{
      residents_xcor(i_ndisp) = lynx_xcor(i), residents_steps(i_ndisp) = lynx_steps(i), residents_ycor(i_ndisp) = lynx_ycor(i), residents_who(i_ndisp) = lynx_who(i), residents_heading(i_ndisp) = lynx_heading(i), residents_prevX(i_ndisp) = lynx_prevX(i), residents_prevY(i_ndisp) = lynx_prevY(i), residents_breed(i_ndisp) = lynx_breed(i), residents_color(i_ndisp) = lynx_color(i), residents_pop(i_ndisp) = lynx_pop(i), residents_sex(i_ndisp) = lynx_sex(i), residents_age(i_ndisp) = lynx_age(i), residents_status(i_ndisp) = lynx_status(i), residents_lastDispX(i_ndisp) = lynx_lastDispX(i), residents_lastDispY(i_ndisp) = lynx_lastDispY(i), residents_nMat(i_ndisp) = lynx_nMat(i), residents_maleID(i_ndisp) = lynx_maleID(i), residents_nFem(i_ndisp) = lynx_nFem(i), residents_rdMortTerr(i_ndisp) = lynx_rdMortTerr(i);
      i_ndisp++;
    }
  }
  //return dispersers_xcor;
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // work on disperser from now on /////////////////////////////////////////////////////////////////
  
  int maxDisp = max(dispersers_steps);
  int step_count = int_0;
  for(int step = 1; step < (maxDisp+1); step++){ //maxDisp
    
    step_count++;
    
    Rcout << "Rcout 1" << std::endl << step_count << std::endl;
    
    // sanity loop
    // if(dispersers_who.size() != nDispLeft){
    //   List L_return = List::create(Named("where") = "beginning_step_2",
    //                                _["Lynx_who"] = lynx_who,
    //                                _["residents_who"] = residents_who,
    //                                _["dispersers_who"] = dispersers_who,
    //                                _["dispersers_steps"] = dispersers_steps,
    //                                _["nLynx"] = nLynx,
    //                                _["nDispLeft"] = nDispLeft,
    //                                _["nRes"] = nRes,
    //                                _["lynx_lastDispX"] = lynx_lastDispX,
    //                                _["deathRoadOld"] = deathRoadOld,
    //                                _["step_count++"] = step_count);
    //   return L_return;
    // }
    
    
    //////////////////////////////////////////////////////////////////////
    /// work on dispersion env around indiv and type of movement
    if(nDispLeft>0){
      //mat with X / Y / id ind / Habitat cell
      IntegerVector CellsDisp_lastDispX(nDispLeft * 9);
      IntegerVector CellsDisp_lastDispY(nDispLeft * 9);
      IntegerVector CellsDisp_pxcor(nDispLeft * 9);
      IntegerVector CellsDisp_pycor(nDispLeft * 9);
      IntegerVector CellsDisp_pxcorHere(nDispLeft * 9);
      IntegerVector CellsDisp_pycorHere(nDispLeft * 9);
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
            CellsDisp_who(i * 9 + l1 + l2 * 3) = dispersers_who(i);
            CellsDisp_steps(i * 9 + l1 + l2 * 3) = dispersers_steps(i);
            CellsDisp_nMat(i * 9 + l1 + l2 * 3) = dispersers_nMat(i);
            CellsDisp_heading(i * 9 + l1 + l2 * 3) = dispersers_heading(i);
            CellsDisp_lastDispX(i * 9 + l1 + l2 * 3) = dispersers_lastDispX(i) - 1 + l1;// because want square around indiv
            CellsDisp_lastDispY(i * 9 + l1 + l2 * 3) = dispersers_lastDispY(i) - 1 + l1;
            CellsDisp_pxcor(i * 9 + l1 + l2 * 3) = dispersers_lastDispX(i) - 1 + l1;// because want square around indiv
            CellsDisp_pycor(i * 9 + l1 + l2 * 3) = dispersers_lastDispY(i) - 1 + l1;
            CellsDisp_pxcorHere(i * 9 + l1 + l2 * 3) = dispersers_lastDispX(i) - 1 + l1;// because want square around indiv
            CellsDisp_pycorHere(i * 9 + l1 + l2 * 3) = dispersers_lastDispY(i) - 1 + l1;// THERE IS A lOT oS SUCCESSIVE RENAMING IN SARAH'S SCRIPT
            CellsDisp_ind(i * 9 + l1 + l2 * 3) = i;
            CellsDisp_hab(i * 9 + l1 + l2 * 3) = HabitatMap(CellsDisp_lastDispY(i * 9 + l1 + l2 * 3), CellsDisp_lastDispX(i * 9 + l1 + l2 * 3));
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
        
        Mat_Chosen(i) = R::rbinom(int_1, pMat * HabFreqDisp_matrix(i)); // number of samples is argument 2, ie matrix type
        if(HabFreqDisp_matrix(i) == 9){ // ie if matrix picked 9 times out of nine, you have a prob one of using matrix type hab
          Mat_Chosen(i) = 1;
        }
        // Sanity check, give resuts just before hend loop
        // if(i == (nDispLeft - 1)){
        //   List L_return = List::create(Named("where") = "hab_freq/mat_chosen",
        //                                _["Lynx_who"] = lynx_who,
        //                                _["step_count++"] = step_count,
        //                                _["CellsDisp_lastDispX"] = CellsDisp_lastDispX,
        //                                _["Mat_Chosen"] = Mat_Chosen);
        //   return L_return;
        // }
        
      }
      
      Rcout << "Rcout 2 matchosen" << std::endl << step_count << std::endl;
      
      //sanity check
      // if(step==5){
      //   List L_return = List::create(Named("where") = "mat_chosen",
      //                                _["step_count++"] = step_count,
      //                                _["dispersers_who"] = dispersers_who,
      //                                _["Mat_Chosen"] = Mat_Chosen);
      //   return L_return;
      // } issue below
      
      // find final number cell for dispersing individuals, taking matrix or disprep only
      int nCellsDispLeft = 0;
      for(int nc = 0; nc<CellsDisp_lastDispX.size(); nc++){ // nc number of cells
        if((Mat_Chosen(CellsDisp_ind(nc)) == 1) & (CellsDisp_hab(nc) == 2)){
          nCellsDispLeft++;
        }
        if((Mat_Chosen(CellsDisp_ind(nc)) == 0) & ((CellsDisp_hab(nc) == 3) | (CellsDisp_hab(nc) == 4))){
          nCellsDispLeft++;
        }
      }
      
      // loop to pick only the right habitat for further dispersion based on whether mat_Chosen = 1
      IntegerVector nextCellsType_lastDispX(nCellsDispLeft);
      IntegerVector nextCellsType_lastDispY(nCellsDispLeft);
      IntegerVector nextCellsType_pxcor(nCellsDispLeft);
      IntegerVector nextCellsType_pycor(nCellsDispLeft);
      IntegerVector nextCellsType_pxcorHere(nCellsDispLeft);
      IntegerVector nextCellsType_pycorHere(nCellsDispLeft);
      IntegerVector nextCellsType_ind(nCellsDispLeft);
      IntegerVector nextCellsType_hab(nCellsDispLeft);
      IntegerVector nextCellsType_who(nCellsDispLeft);
      IntegerVector nextCellsType_steps(nCellsDispLeft);
      IntegerVector nextCellsType_nMat(nCellsDispLeft);
      IntegerVector nextCellsType_heading(nCellsDispLeft);
      for(int i = 0, p = 0; i<CellsDisp_ind.size(); i++){
        if((Mat_Chosen(CellsDisp_ind(i)) == 1) & (CellsDisp_hab(i) == 2)){
          //stop("Went into loop");
          nextCellsType_lastDispX(p) = CellsDisp_lastDispX(i);
          nextCellsType_lastDispY(p) = CellsDisp_lastDispY(i);
          nextCellsType_pxcor(p) = CellsDisp_pxcor(i);
          nextCellsType_pycor(p) = CellsDisp_pycor(i);
          nextCellsType_pxcorHere(p) = CellsDisp_pxcorHere(i);
          nextCellsType_pycorHere(p) = CellsDisp_pycorHere(i);
          nextCellsType_ind(p) = CellsDisp_ind(i);
          nextCellsType_hab(p) = CellsDisp_hab(i);
          nextCellsType_who(p) = CellsDisp_who(i);
          nextCellsType_steps(p) = CellsDisp_steps(i);
          nextCellsType_nMat(p) = CellsDisp_nMat(i);
          nextCellsType_heading(p) = CellsDisp_heading(i);
          p++;
        }
        if((Mat_Chosen(CellsDisp_ind(i)) == 0) & ((CellsDisp_hab(i) == 3) | (CellsDisp_hab(i) == 4))){
          nextCellsType_lastDispX(p) = CellsDisp_lastDispX(i);
          nextCellsType_lastDispY(p) = CellsDisp_lastDispY(i);
          nextCellsType_pxcor(p) = CellsDisp_pxcor(i);
          nextCellsType_pycor(p) = CellsDisp_pycor(i);
          nextCellsType_pxcorHere(p) = CellsDisp_pxcorHere(i);
          nextCellsType_pycorHere(p) = CellsDisp_pycorHere(i);
          nextCellsType_ind(p) = CellsDisp_ind(i);
          nextCellsType_hab(p) = CellsDisp_hab(i);
          nextCellsType_who(p) = CellsDisp_who(i);
          nextCellsType_steps(p) = CellsDisp_steps(i);
          nextCellsType_nMat(p) = CellsDisp_nMat(i);
          nextCellsType_heading(p) = CellsDisp_heading(i);
          p++;
        }
      }// end second ind loop
      
      Rcout << "Rcout 3 nextCellsType" << std::endl << step_count << std::endl;
      
      // final steps elements def here for push back
      IntegerVector ChosenCells_lastDispX(0);
      IntegerVector ChosenCells_lastDispY(0);
      IntegerVector ChosenCells_pxcor(0);
      IntegerVector ChosenCells_pycor(0);
      IntegerVector ChosenCells_pxcorHere(0);
      IntegerVector ChosenCells_pycorHere(0);
      IntegerVector ChosenCells_hab(0);
      IntegerVector ChosenCells_ind(0);
      IntegerVector ChosenCells_who(0);
      IntegerVector ChosenCells_steps(0);
      IntegerVector ChosenCells_IsMoveCorr(0);
      IntegerVector ChosenCells_nMat(0);
      IntegerVector ChosenCells_heading(0);
      
      // part on potential correlation in movements, none on first move but then some/////////////////////////////////////////
      
      //step++;
      IntegerVector randLines_move = sample((nextCellsType_hab.size()), (nextCellsType_hab.size()), false) - 1; // number of samples is argument 2;
      if(step == 1){
        // here randomly shuffling line before picking one per ind, done above
        for(int ind = 0; ind < nDispLeft; ind++){
          double p = 0.5;
          while(p<1){
            for(int l = 0; l<randLines_move.size(); l++){
              //if(nextCellsType_ind(randLines_move(l) == ind)){ // here keeps last one, while was not working
              if((nextCellsType_ind(randLines_move(l)) == ind) & (p<1)){
                ChosenCells_lastDispX.push_back(nextCellsType_lastDispX(randLines_move(l)));
                ChosenCells_lastDispY.push_back(nextCellsType_lastDispY(randLines_move(l)));
                ChosenCells_pxcor.push_back(nextCellsType_pxcor(randLines_move(l)));
                ChosenCells_pycor.push_back(nextCellsType_pycor(randLines_move(l)));
                ChosenCells_pxcorHere.push_back(nextCellsType_pxcorHere(randLines_move(l)));
                ChosenCells_pycorHere.push_back(nextCellsType_pycorHere(randLines_move(l)));
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
        // sanity check
        // List L_return = List::create(Named("nCellsDispLeft") = nCellsDispLeft,
        //                               _["nDispLeft"] = nDispLeft,
        //                               //_["randLines_move"] = randLines_move,
        //                               _["ChosenCell_lastDispX"] = ChosenCells_lastDispX,
        //                               _["ChosenCell_who"] = ChosenCells_who,
        //                               _["ChosenCell_ind"] = ChosenCells_ind);
        //  return  L_return;
        
      } //Works up to here
      
      Rcout << "Rcout 4 after steps 0" << std::endl << step_count << std::endl;
      
      // sanity check
      // List L_return = List::create(Named("where") = "nCellsDispLeft done",
      // _["nCellsDispLeft"] = nCellsDispLeft
      //                              _["nDispLeft"] = nDispLeft,
      //                              _["randLines_move"] = randLines_move,
      //                              _["ChosenCell_lastDispX"] = ChosenCells_lastDispX,
      //                              _["ChosenCell_who"] = ChosenCells_who);
      // return  L_return;
      
      if(step>1){ // i.e if step more than first //////////////////////////////////////////////////
        
        //Have to define here and push_back values to be able to use latter on outside loop where it is filled up
        // final move of indiv with no correlated movements
        IntegerVector ChosenCellsNoCorr_who(0);
        IntegerVector ChosenCellsNoCorr_ind(0);
        IntegerVector ChosenCellsNoCorr_hab(0);
        IntegerVector ChosenCellsNoCorr_lastDispX(0);
        IntegerVector ChosenCellsNoCorr_lastDispY(0);
        IntegerVector ChosenCellsNoCorr_pxcor(0);
        IntegerVector ChosenCellsNoCorr_pycor(0);
        IntegerVector ChosenCellsNoCorr_pxcorHere(0);
        IntegerVector ChosenCellsNoCorr_pycorHere(0);
        IntegerVector ChosenCellsNoCorr_steps(0);
        IntegerVector ChosenCellsNoCorr_nMat(0);
        IntegerVector ChosenCellsNoCorr_heading(0);
        IntegerVector ChosenCellsNoCorr_IsMoveCorr(0);
        // final move indiv with correlated movement
        IntegerVector ChosenCellsYesCorr_who(0);
        IntegerVector ChosenCellsYesCorr_ind(0);
        IntegerVector ChosenCellsYesCorr_hab(0);
        IntegerVector ChosenCellsYesCorr_prefDir(0);
        IntegerVector ChosenCellsYesCorr_pxcor(0);
        IntegerVector ChosenCellsYesCorr_pycor(0);
        IntegerVector ChosenCellsYesCorr_pxcorHere(0);
        IntegerVector ChosenCellsYesCorr_pycorHere(0);
        IntegerVector ChosenCellsYesCorr_lastDispX(0);
        IntegerVector ChosenCellsYesCorr_lastDispY(0);
        IntegerVector ChosenCellsYesCorr_steps(0);
        IntegerVector ChosenCellsYesCorr_nMat(0);
        IntegerVector ChosenCellsYesCorr_heading(0);
        IntegerVector ChosenCellsYesCorr_IsMoveCorr(0);
        
        // sort by ind and add current position next to potential ones in nextCellsType vectors
        IntegerVector nextCellsType_indSortIndex = IntOrderIndex(nextCellsType_ind);
        IntegerVector nextCellsType_indSorted = nextCellsType_ind.size();
        IntegerVector nextCellsType_habSorted = nextCellsType_ind.size();
        IntegerVector nextCellsType_pxcorSorted = nextCellsType_ind.size();
        IntegerVector nextCellsType_pycorSorted = nextCellsType_ind.size();
        IntegerVector nextCellsType_pxcorHereSorted = nextCellsType_ind.size();
        IntegerVector nextCellsType_pycorHereSorted = nextCellsType_ind.size();
        IntegerVector nextCellsType_lastDispXSorted = nextCellsType_ind.size();
        IntegerVector nextCellsType_lastDispYSorted = nextCellsType_ind.size();
        IntegerVector nextCellsType_whoSorted = nextCellsType_who.size();
        IntegerVector nextCellsType_stepsSorted = nextCellsType_steps.size();
        IntegerVector nextCellsType_nMatSorted = nextCellsType_nMat.size();
        IntegerVector nextCellsType_headingSorted = nextCellsType_heading.size();
        
        for(int l = 0; l< nextCellsType_ind.size(); l++){
          nextCellsType_indSorted(l) = nextCellsType_ind(nextCellsType_indSortIndex(l));
          nextCellsType_pxcorSorted(l) = nextCellsType_pxcor(nextCellsType_indSortIndex(l));
          nextCellsType_pycorSorted(l) = nextCellsType_pycor(nextCellsType_indSortIndex(l));
          nextCellsType_lastDispXSorted(l) = nextCellsType_lastDispX(nextCellsType_indSortIndex(l));
          nextCellsType_lastDispYSorted(l) = nextCellsType_lastDispY(nextCellsType_indSortIndex(l));
          nextCellsType_habSorted(l) = nextCellsType_hab(nextCellsType_indSortIndex(l));
          nextCellsType_whoSorted(l) = nextCellsType_who(nextCellsType_indSortIndex(l));
          nextCellsType_stepsSorted(l) = nextCellsType_steps(nextCellsType_indSortIndex(l));
          nextCellsType_nMatSorted(l) = nextCellsType_nMat(nextCellsType_indSortIndex(l));
          nextCellsType_headingSorted(l) = nextCellsType_heading(nextCellsType_indSortIndex(l));
          for(int j = 0; j<dispersers_who.size(); j++){
            if(dispersers_who(j) == nextCellsType_whoSorted(l)){
              nextCellsType_pxcorHereSorted(l) = dispersers_lastDispX(j);
              nextCellsType_pycorHereSorted(l) = dispersers_lastDispY(j);
            }
          }
        }
        
        Rcout << "Rcout 5  nextCellsType step>0" << std::endl << step_count << std::endl;
        
        // dispersers with some steps left, ie final matrix
        IntegerVector WStepsLeft = WhichAbove(nextCellsType_stepsSorted, -10); //before was step - 1 but now all can go, selected before//- 1 because I want >= behavior from > function
        IntegerVector nextCellsType_indF = IntVecSubIndex(nextCellsType_indSorted, WStepsLeft);
        IntegerVector nextCellsType_habF = IntVecSubIndex(nextCellsType_habSorted, WStepsLeft);
        IntegerVector nextCellsType_pxcorF = IntVecSubIndex(nextCellsType_pxcorSorted, WStepsLeft);
        IntegerVector nextCellsType_pycorF = IntVecSubIndex(nextCellsType_pycorSorted, WStepsLeft);
        IntegerVector nextCellsType_pxcorHereF = IntVecSubIndex(nextCellsType_pxcorHereSorted, WStepsLeft);
        IntegerVector nextCellsType_pycorHereF = IntVecSubIndex(nextCellsType_pycorHereSorted, WStepsLeft);
        IntegerVector nextCellsType_lastDispXF = IntVecSubIndex(nextCellsType_lastDispXSorted, WStepsLeft);
        IntegerVector nextCellsType_lastDispYF = IntVecSubIndex(nextCellsType_lastDispYSorted, WStepsLeft);
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
          IntegerVector nextCellsTypeNoCorr_habF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_lastDispXF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_lastDispYF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_pxcorF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_pycorF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_pxcorHereF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_pycorHereF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_whoF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_stepsF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_nMatF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_headingF(nCorr0);
          IntegerVector nextCellsTypeNoCorr_IsMoveCorrF(nCorr0);
          for(int l = 0; l<nCorr0;l++){
            nextCellsTypeNoCorr_indF(l) = nextCellsType_indF(noCorr_Lind(l));
            nextCellsTypeNoCorr_pxcorF(l) = nextCellsType_pxcorF(noCorr_Lind(l));
            nextCellsTypeNoCorr_pycorF(l) = nextCellsType_pycorF(noCorr_Lind(l));
            nextCellsTypeNoCorr_pxcorHereF(l) = nextCellsType_pxcorHereF(noCorr_Lind(l));
            nextCellsTypeNoCorr_pycorHereF(l) = nextCellsType_pycorHereF(noCorr_Lind(l));
            nextCellsTypeNoCorr_lastDispXF(l) = nextCellsType_lastDispXF(noCorr_Lind(l));
            nextCellsTypeNoCorr_lastDispYF(l) = nextCellsType_lastDispYF(noCorr_Lind(l));
            nextCellsTypeNoCorr_habF(l) = nextCellsType_habF(noCorr_Lind(l));
            nextCellsTypeNoCorr_whoF(l) = nextCellsType_whoF(noCorr_Lind(l));
            nextCellsTypeNoCorr_stepsF(l) = nextCellsType_stepsF(noCorr_Lind(l));
            nextCellsTypeNoCorr_nMatF(l) = nextCellsType_nMatF(noCorr_Lind(l));
            nextCellsTypeNoCorr_headingF(l) = nextCellsType_headingF(noCorr_Lind(l));
            nextCellsTypeNoCorr_IsMoveCorrF(l) = nextCellsType_IsMoveCorrF(noCorr_Lind(l));
          }
          // now pick just one cell to move to for each ind
          IntegerVector UniqueLinesIndNoCorr = IntPosOneOfEach(nextCellsTypeNoCorr_indF);
          for(int i = 0; i<UniqueLinesIndNoCorr.size(); i++){
            ChosenCellsNoCorr_ind.push_back(nextCellsTypeNoCorr_indF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_pxcor.push_back(nextCellsTypeNoCorr_pxcorF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_pycor.push_back(nextCellsTypeNoCorr_pycorF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_pxcorHere.push_back(nextCellsTypeNoCorr_pxcorHereF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_pycorHere.push_back(nextCellsTypeNoCorr_pycorHereF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_lastDispX.push_back(nextCellsTypeNoCorr_lastDispXF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_lastDispY.push_back(nextCellsTypeNoCorr_lastDispYF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_hab.push_back(nextCellsTypeNoCorr_habF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_who.push_back(nextCellsTypeNoCorr_whoF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_steps.push_back(nextCellsTypeNoCorr_stepsF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_nMat.push_back(nextCellsTypeNoCorr_nMatF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_heading.push_back(nextCellsTypeNoCorr_headingF(UniqueLinesIndNoCorr(i)));
            ChosenCellsNoCorr_IsMoveCorr.push_back(nextCellsTypeNoCorr_IsMoveCorrF(UniqueLinesIndNoCorr(i)));
          }
          // sanity check
          // List L_return = List::create(Named("where") = "ChosenCellsNoCorr",
          //                              _["ChosenCellsNoCorr_who"] = ChosenCellsNoCorr_who,
          //                              _["ChosenCellsNoCorr_ind"] = ChosenCellsNoCorr_ind,
          //                              _["ChosenCellsNoCorr_pycorHere"] = ChosenCellsNoCorr_pycorHere);
          // return  L_return;
        }// end noCorr move indiv, works up to here
        
        Rcout << "Rcout 6 uncor mov" << std::endl << step_count << std::endl;
        
        // now deal with indiv with correlated movement
        if(nCorr1 > int_0){
          //stop("got into loop");
          IntegerVector YesCorr_Lind = WhichEqual(nextCellsType_IsMoveCorrF, int_1);
          //define vector to fill
          IntegerVector nextCellsTypeYesCorr_indF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_pxcorHereF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_pycorHereF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_pxcorF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_pycorF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_lastDispXF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_lastDispYF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_habF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_whoF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_stepsF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_nMatF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_headingF(YesCorr_Lind.size());
          IntegerVector nextCellsTypeYesCorr_IsMoveCorrF(YesCorr_Lind.size());
          for(int l = 0; l<YesCorr_Lind.size(); l++){
            nextCellsTypeYesCorr_indF(l) = nextCellsType_indF(YesCorr_Lind(l));
            nextCellsTypeYesCorr_habF(l) = nextCellsType_habF(YesCorr_Lind(l));
            nextCellsTypeYesCorr_pxcorHereF(l) = nextCellsType_pxcorHereF(YesCorr_Lind(l));
            nextCellsTypeYesCorr_pycorHereF(l) = nextCellsType_pycorHereF(YesCorr_Lind(l));
            nextCellsTypeYesCorr_pxcorF(l) = nextCellsType_pxcorF(YesCorr_Lind(l));
            nextCellsTypeYesCorr_pycorF(l) = nextCellsType_pycorF(YesCorr_Lind(l));
            nextCellsTypeYesCorr_lastDispXF(l) = nextCellsType_lastDispXF(YesCorr_Lind(l));
            nextCellsTypeYesCorr_lastDispYF(l) = nextCellsType_lastDispYF(YesCorr_Lind(l));
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
            nextCellsTypeYesCorr_DirF(l) = towards_simple_unique(nextCellsTypeYesCorr_pxcorHereF(l), nextCellsTypeYesCorr_pycorHereF(l),
                                      nextCellsTypeYesCorr_pxcorF(l), nextCellsTypeYesCorr_pycorF(l));
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
            if((nextCellsTypeYesCorr_pxcorHereF(l) == nextCellsTypeYesCorr_pxcorF(l)) & (nextCellsTypeYesCorr_pycorHereF(l) == nextCellsTypeYesCorr_pycorF(l))){
              nextCellsTypeYesCorr_prefDirF(l) = int_3;
            }
          }
          
          Rcout << "Rcout 7 cor mov" << std::endl << step_count << std::endl;
          
          // keep only one line per indiv, with lowest value of rank
          IntegerVector unique_nextCellsTypeYesCorr_whoF = unique(nextCellsTypeYesCorr_whoF);
          unique_nextCellsTypeYesCorr_whoF = sortInt(unique_nextCellsTypeYesCorr_whoF);
          // subset for one of the lower prefdir values
          for(int i = 0; i<unique_nextCellsTypeYesCorr_whoF.size(); i++){
            ChosenCellsYesCorr_prefDir.push_back(int_100);// set to 100 to be able to replace by values within table
            ChosenCellsYesCorr_who.push_back(int_100);
            ChosenCellsYesCorr_ind.push_back(int_100);
            ChosenCellsYesCorr_hab.push_back(int_100);
            ChosenCellsYesCorr_pxcorHere.push_back(int_100);
            ChosenCellsYesCorr_pycorHere.push_back(int_100);
            ChosenCellsYesCorr_pxcor.push_back(int_100);
            ChosenCellsYesCorr_pycor.push_back(int_100);
            ChosenCellsYesCorr_lastDispX.push_back(int_100);
            ChosenCellsYesCorr_lastDispY.push_back(int_100);
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
                ChosenCellsYesCorr_hab(i) = nextCellsTypeYesCorr_habF(l);
                ChosenCellsYesCorr_pxcorHere(i) = nextCellsTypeYesCorr_pxcorHereF(l);
                ChosenCellsYesCorr_pycorHere(i) = nextCellsTypeYesCorr_pycorHereF(l);
                ChosenCellsYesCorr_lastDispX(i) = nextCellsTypeYesCorr_lastDispXF(l);
                ChosenCellsYesCorr_lastDispY(i) = nextCellsTypeYesCorr_lastDispYF(l);
                ChosenCellsYesCorr_pxcor(i) = nextCellsTypeYesCorr_pxcorF(l);
                ChosenCellsYesCorr_pycor(i) = nextCellsTypeYesCorr_pycorF(l);
                ChosenCellsYesCorr_steps(i) = nextCellsTypeYesCorr_stepsF(l);
                ChosenCellsYesCorr_nMat(i) = nextCellsTypeYesCorr_nMatF(l);
                ChosenCellsYesCorr_heading(i) = nextCellsTypeYesCorr_headingF(l);
                ChosenCellsYesCorr_IsMoveCorr(i) = nextCellsTypeYesCorr_IsMoveCorrF(l);
              }
            }
          }
        }
        
        // sanity check
        // List L_return = List::create(Named("where") = "correlated move part",
        //                              _["nextCellsType_IsMoveCorrF"] = nextCellsType_IsMoveCorrF,
        //                              _["nCorr1"] = nCorr1,
        //                              _["ChosenCellsYesCorr_prefDir"] = ChosenCellsYesCorr_prefDir,
        //                              _["ChosenCellsYesCorr_who"] = ChosenCellsYesCorr_who,
        //                              _["ChosenCellsNoCorr_who"] = ChosenCellsNoCorr_who,
        //                              _["ChosenCellsYesCorr_hab"] = ChosenCellsYesCorr_hab);
        // return L_return;
        
        // get together moves for correlated and uncorrelated individuals
        int nLChosenCellsYesCorr = ChosenCellsYesCorr_who.size();
        int nLChosenCellsNoCorr = ChosenCellsNoCorr_who.size(); //ChosenCellsNoCorr
        int nChosenCells = nLChosenCellsYesCorr + nLChosenCellsNoCorr;
        
        //sanity check
        //if(step == 10){
        // if(nDispLeft > (nLChosenCellsYesCorr + nLChosenCellsNoCorr)){
        //   List L_return = List::create(Named("where") = "Yes and No corr move",
        //                                _["step"] = step,
        //                                _["nCellsDispLeft"] = nCellsDispLeft,
        //                                _["nDispLeftOld"] = nDispLeftOld,
        //                                _["nDispLeft"] = nDispLeft,
        //                                _["nextCellsType_IsMoveCorrF"] = nextCellsType_IsMoveCorrF,
        //                                _["ChosenCellsYesCorr_who"] = ChosenCellsYesCorr_who,
        //                                _["ChosenCellsNoCorr_who"] = ChosenCellsNoCorr_who,
        //                                _["CellsDisp_who"] = CellsDisp_who,
        //                                _["nCellsDispLeft"] = nCellsDispLeft,
        //                                _["nextCellsType_who"] = nextCellsType_who,
        //                                _["nextCellsType_whoSorted"] = nextCellsType_whoSorted,
        //                                _["nLChosenCellsYesCorr"] = nLChosenCellsYesCorr,
        //                                _["nLChosenCellsNoCorr"] = nLChosenCellsNoCorr,
        //                                _["Dispersers_who"] = dispersers_who,
        //                                _["ChosenCell_who"] = ChosenCells_who,
        //                                _["Dispersers_steps"] = ChosenCells_steps,
        //                                _["ChosenCell_steps"] = ChosenCells_steps,
        //                                _["deathRoadOld"] = deathRoadOld);
        //   return  L_return;
        // }
        
        Rcout << "Rcout 8 cor mov2" << std::endl << step_count << std::endl;
        
        for(int l = 0; l<nChosenCells; l++){
          if(l < nLChosenCellsYesCorr){
            ChosenCells_who.push_back(ChosenCellsYesCorr_who(l));
            ChosenCells_ind.push_back(ChosenCellsYesCorr_ind(l));
            ChosenCells_hab.push_back(ChosenCellsYesCorr_hab(l));
            ChosenCells_pxcorHere.push_back(ChosenCellsYesCorr_pxcorHere(l));
            ChosenCells_pycorHere.push_back(ChosenCellsYesCorr_pycorHere(l));
            ChosenCells_pxcor.push_back(ChosenCellsYesCorr_pxcor(l));
            ChosenCells_pycor.push_back(ChosenCellsYesCorr_pycor(l));
            ChosenCells_lastDispX.push_back(ChosenCellsYesCorr_lastDispX(l));
            ChosenCells_lastDispY.push_back(ChosenCellsYesCorr_lastDispY(l));
            ChosenCells_nMat.push_back(ChosenCellsYesCorr_nMat(l));
            ChosenCells_heading.push_back(ChosenCellsYesCorr_heading(l));
            ChosenCells_IsMoveCorr.push_back(ChosenCellsYesCorr_IsMoveCorr(l));
            ChosenCells_steps.push_back(ChosenCellsYesCorr_steps(l));
          }else{
            ChosenCells_who.push_back(ChosenCellsNoCorr_who(l - nLChosenCellsYesCorr));
            ChosenCells_ind.push_back(ChosenCellsNoCorr_ind(l - nLChosenCellsYesCorr));
            ChosenCells_hab.push_back(ChosenCellsNoCorr_hab(l - nLChosenCellsYesCorr));
            ChosenCells_pxcorHere.push_back(ChosenCellsNoCorr_pxcorHere(l - nLChosenCellsYesCorr));
            ChosenCells_pycorHere.push_back(ChosenCellsNoCorr_pycorHere(l - nLChosenCellsYesCorr));
            ChosenCells_pxcor.push_back(ChosenCellsNoCorr_pxcor(l - nLChosenCellsYesCorr));
            ChosenCells_pycor.push_back(ChosenCellsNoCorr_pycor(l - nLChosenCellsYesCorr));
            ChosenCells_lastDispX.push_back(ChosenCellsNoCorr_lastDispX(l - nLChosenCellsYesCorr));
            ChosenCells_lastDispY.push_back(ChosenCellsNoCorr_lastDispY(l - nLChosenCellsYesCorr));
            ChosenCells_nMat.push_back(ChosenCellsNoCorr_nMat(l - nLChosenCellsYesCorr));
            ChosenCells_heading.push_back(ChosenCellsNoCorr_heading(l - nLChosenCellsYesCorr));
            ChosenCells_IsMoveCorr.push_back(ChosenCellsNoCorr_IsMoveCorr(l - nLChosenCellsYesCorr));
            ChosenCells_steps.push_back(ChosenCellsNoCorr_steps(l - nLChosenCellsYesCorr));
          }
        }
      }
      
      //sanity
      //if(step == 10){
      // if(nDispLeft > ChosenCells_who.size()){
      //   List L_return = List::create(Named("where") = "chosenCells done",
      //                                _["step"] = step,
      //                                _["nCellsDispLeft"] = nCellsDispLeft,
      //                                _["nDispLeftOld"] = nDispLeftOld,
      //                                _["nDispLeft"] = nDispLeft,
      //                                _["Dispersers_who"] = dispersers_who,
      //                                _["ChosenCell_lastDispX"] = ChosenCells_lastDispX,
      //                                _["ChosenCell_who"] = ChosenCells_who,
      //                                _["ChosenCell_steps"] = ChosenCells_steps,
      //                                _["deathRoadOld"] = deathRoadOld);
      //   return  L_return;
      // }
      
      // work on chosenMat and chosenDisp matrices and their processing /////////////////////////////////////////
      IntegerVector MatInd = WhichEqual(ChosenCells_hab, int_2);
      int nMatInd = MatInd.size();
      int nDispInd = ChosenCells_nMat.size() - nMatInd;
      IntegerVector ChosenMat_who(nMatInd);
      IntegerVector ChosenMat_ind(nMatInd);
      IntegerVector ChosenMat_hab(nMatInd);
      IntegerVector ChosenMat_pxcorHere(nMatInd);
      IntegerVector ChosenMat_pycorHere(nMatInd);
      IntegerVector ChosenMat_pxcor(nMatInd);
      IntegerVector ChosenMat_pycor(nMatInd);
      IntegerVector ChosenMat_lastDispX(nMatInd);
      IntegerVector ChosenMat_lastDispY(nMatInd);
      IntegerVector ChosenMat_nMat(nMatInd);
      IntegerVector ChosenMat_heading(nMatInd);
      // disp
      IntegerVector ChosenDisp_who(nDispInd);
      IntegerVector ChosenDisp_ind(nDispInd);
      IntegerVector ChosenDisp_hab(nDispInd);
      IntegerVector ChosenDisp_pxcorHere(nDispInd);
      IntegerVector ChosenDisp_pycorHere(nDispInd);
      IntegerVector ChosenDisp_pxcor(nDispInd);
      IntegerVector ChosenDisp_pycor(nDispInd);
      IntegerVector ChosenDisp_lastDispX(nDispInd);
      IntegerVector ChosenDisp_lastDispY(nDispInd);
      IntegerVector ChosenDisp_nMat(nDispInd);
      IntegerVector ChosenDisp_heading(nDispInd);
      for(int l = 0, p = 0, q = 0; l<ChosenCells_nMat.size(); l++){
        if(ChosenCells_hab(l) == 2){
          if(nMatInd>0){
            ChosenMat_who(p) = ChosenCells_who(l);
            ChosenMat_ind(p) = ChosenCells_ind(l);
            ChosenMat_hab(p) = ChosenCells_hab(l);
            ChosenMat_pxcorHere(p) = ChosenCells_pxcorHere(l);
            ChosenMat_pycorHere(p) = ChosenCells_pycorHere(l);
            ChosenMat_pxcor(p) = ChosenCells_pxcor(l);
            ChosenMat_pycor(p) = ChosenCells_pycor(l);
            ChosenMat_lastDispX(p) = ChosenCells_lastDispX(l);
            ChosenMat_lastDispY(p) = ChosenCells_lastDispY(l);
            ChosenMat_nMat(p) = ChosenCells_nMat(l) + 1;
            // code memory in movement bit here
            if((ChosenMat_nMat(p) + 1) == nMatMax){
              ChosenMat_pxcor(p) = ChosenMat_pxcor(p); // there is really that
              ChosenMat_pycor(p) = ChosenMat_pycor(p);
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
            ChosenDisp_pxcor(q) = ChosenCells_pxcor(l);
            ChosenDisp_pycor(q) = ChosenCells_pycor(l);
            ChosenDisp_pxcorHere(q) = ChosenCells_pxcorHere(l);
            ChosenDisp_pycorHere(q) = ChosenCells_pycorHere(l);
            ChosenDisp_lastDispX(q) = ChosenCells_lastDispX(l);
            ChosenDisp_lastDispY(q) = ChosenCells_lastDispY(l);
            ChosenDisp_nMat(q) = 0;
            q++;
          }
        }
      }
      
      Rcout << "Rcout 9 chosenDisp" << std::endl << step_count << std::endl;
      
      // process ChosenMat lynxMemory already integrated within loops above
      // same for pxcor /lastdisp update
      
      // now reupdate chosen cells matrix by binding ChosenDisp and ChosenMat
      for(int l = 0; l<ChosenCells_who.size(); l++){
        if(l<nDispInd){
          ChosenCells_who(l) = ChosenDisp_who(l);
          ChosenCells_ind(l) = ChosenDisp_ind(l);
          ChosenCells_hab(l) =ChosenDisp_hab(l);
          ChosenCells_pxcorHere(l) = ChosenDisp_pxcorHere(l);
          ChosenCells_pycorHere(l) = ChosenDisp_pycorHere(l);
          ChosenCells_pxcor(l) = ChosenDisp_pxcor(l);
          ChosenCells_pycor(l) = ChosenDisp_pycor(l);
          ChosenCells_lastDispX(l) = ChosenDisp_lastDispX(l);
          ChosenCells_lastDispY(l) = ChosenDisp_lastDispY(l);
          ChosenCells_nMat(l) = ChosenDisp_nMat(l);
        }
        else{
          ChosenCells_who(l) = ChosenMat_who(l - nDispInd);
          ChosenCells_ind(l) = ChosenMat_ind(l - nDispInd);
          ChosenCells_hab(l) = ChosenMat_hab(l - nDispInd);
          ChosenCells_pxcorHere(l) = ChosenMat_pxcorHere(l - nDispInd);
          ChosenCells_pycorHere(l) = ChosenMat_pycorHere(l - nDispInd);
          ChosenCells_pxcor(l) = ChosenMat_pxcor(l - nDispInd);
          ChosenCells_pycor(l) = ChosenMat_pycor(l - nDispInd);
          ChosenCells_lastDispX(l) = ChosenMat_lastDispX(l - nDispInd);
          ChosenCells_lastDispY(l) = ChosenMat_lastDispY(l - nDispInd);
          ChosenCells_nMat(l) = ChosenMat_nMat(l - nDispInd);
        }
      }
      //update connectivity map with +1 when dipserser stp on cell
      for(int l = 0; l<ChosenCells_pxcor.size(); l++){
        connectivityMap(ChosenCells_pycor(l) , ChosenCells_pxcor(l)) +=1;
      }
      
      Rcout << "Rcout 10 chosenCells" << std::endl << step_count << std::endl;
      
      //////////////////////////////////////////////////////////////////////////////////////////////
      // bit on road mortality
      IntegerVector deathRoad(ChosenCells_pxcor.size());
      for(int l = 0; l<ChosenCells_pxcor.size(); l++){
        deathRoad(l) = R::rbinom(int_1, (roadMortMap(ChosenCells_pycor(l), ChosenCells_pxcor(l)) / corrFactorDisp));
        if(floorTimeSim == startSimYear){ // cannot die first year
          deathRoad(l) = int_0;
        }
      }
      ncoll_ncoll.push_back(sum(deathRoad));
      ncoll_time.push_back(floorTimeSim);
      for(int l = 0; l<ChosenCells_pxcor.size(); l++){
        if(roadMortMap(ChosenCells_pycor(l), ChosenCells_pxcor(l)) == 1){ // force death on border to simulate emigration but does not count as death in line above
          deathRoad(l) = int_1;
        }
      }
      //save some data on dead individual
      for(int l = 0; l < deathRoad.size(); l++){
        if(deathRoad(l) == int_1){
          deadLynxColl.push_back(ChosenCells_who(l), "who");
          deadLynxColl.push_back(ChosenCells_nMat(l), "heading");
          //deadLynxColl.push_back(ChosenCells_who(l), "who"); unusefull duplicate
          deadLynxColl.push_back(ChosenCells_steps(l), "steps");
          deadLynxColl.push_back(ChosenCells_lastDispX(l), "lastDispX");
          deadLynxColl.push_back(ChosenCells_lastDispY(l), "lastDispY");
        }
      }
      // complete deadDisp
      deathRoadOld = clone(deathRoad); // to track issue in next loop
      int deadDispLine = WhichEqual(deadDisp["time"], floorTimeSim)(0);
      IntegerVector deadDisp_nDispDeadColl = deadDisp["nDispDeadColl"];
      deadDisp_nDispDeadColl(deadDispLine) = sum(deathRoad);
      deadDisp["nDispDeadColl"] = deadDisp_nDispDeadColl;
      
      Rcout << "Rcout 11 deadDisp" << std::endl << step_count << std::endl;
      
      // create dispdispersers_who_newerser new table, first two lines were some test
      //deathRoad(0) = 1, deathRoad(1) = 1, deathRoad(2) = 1, deathRoad(3) = 1 , deathRoad(4) = 1;
      //IntegerVector rand_values = sample(5, 10, true) - 1; // number of samples is argument 2; -1 to make it start at 0
      IntegerVector who_ord = IntOrderIndex(ChosenCells_who);
      IntegerVector alive_who_ord = WhichEqual(deathRoad[who_ord], int_0);
      IntegerVector pos_alive_who_ord = who_ord[alive_who_ord];
      // values to pick from ChosenCells
      IntegerVector dispersers_who_new = ChosenCells_who[pos_alive_who_ord], dispersers_ind_new = ChosenCells_ind[pos_alive_who_ord], dispersers_hab_new = ChosenCells_hab[pos_alive_who_ord], dispersers_pxcorHere_new = ChosenCells_pxcorHere[pos_alive_who_ord], dispersers_pycorHere_new = ChosenCells_pycorHere[pos_alive_who_ord], dispersers_pxcor_new = ChosenCells_pxcor[pos_alive_who_ord], dispersers_pycor_new = ChosenCells_pycor[pos_alive_who_ord], dispersers_lastDispX_new = ChosenCells_lastDispX[pos_alive_who_ord], dispersers_lastDispY_new = ChosenCells_lastDispY[pos_alive_who_ord], dispersers_nMat_new = ChosenCells_nMat[pos_alive_who_ord];
      IntegerVector dispersers_steps_new = ChosenCells_steps[pos_alive_who_ord]; // added here because was missing for searchterritory
      int nDisp_new = dispersers_who_new.size();
      // Values to get back from original dispersers data
      IntegerVector index_dispersers_dispersers_new(dispersers_who_new.size());
      if(nDisp_new >= int_1){
        for(int i = 0; i<nDisp_new; i++){
          double p = 0.5;
          while(p<1){
            for(int j = 0; j<dispersers_who.size(); j++){
              if(dispersers_who(j) == dispersers_who_new(i) )
                index_dispersers_dispersers_new(i) = j;
              p += 1;
            }
          }
        }
      }
      IntegerVector dispersers_xcor_new = dispersers_xcor[index_dispersers_dispersers_new], dispersers_ycor_new = dispersers_ycor[index_dispersers_dispersers_new], dispersers_heading_new = dispersers_heading[index_dispersers_dispersers_new], dispersers_prevX_new = dispersers_prevX[index_dispersers_dispersers_new], dispersers_prevY_new = dispersers_prevY[index_dispersers_dispersers_new], dispersers_age_new = dispersers_age[index_dispersers_dispersers_new], dispersers_maleID_new = dispersers_maleID[index_dispersers_dispersers_new], dispersers_nFem_new = dispersers_nFem[index_dispersers_dispersers_new], dispersers_rdMortTerr_new = dispersers_rdMortTerr[index_dispersers_dispersers_new];
      StringVector dispersers_breed_new = dispersers_breed[index_dispersers_dispersers_new], dispersers_color_new = dispersers_color[index_dispersers_dispersers_new], dispersers_pop_new = dispersers_pop[index_dispersers_dispersers_new], dispersers_sex_new = dispersers_sex[index_dispersers_dispersers_new], dispersers_status_new = dispersers_status[index_dispersers_dispersers_new];
      
      Rcout << "Rcout 12 dispersers new" << std::endl << step_count << std::endl;
      
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // update all lynx values //////////////////////////////////////////////////////
      
      // get size of new table, made of residents and dispersers still alive, codes as _new
      int nLynx_new = nDisp_new + nRes;
      IntegerVector res_index(nRes);
      for(int i = 0; i<nRes; i++){
        res_index(i) = i;
      }
      // Initiate the vector for new lynx data
      IntegerVector lynx_xcor_new(nLynx_new), lynx_steps_new(nLynx_new), lynx_ycor_new(nLynx_new), lynx_who_new(nLynx_new), lynx_heading_new(nLynx_new), lynx_prevX_new(nLynx_new), lynx_prevY_new(nLynx_new), lynx_age_new(nLynx_new), lynx_lastDispX_new(nLynx_new), lynx_lastDispY_new(nLynx_new), lynx_nMat_new(nLynx_new), lynx_maleID_new(nLynx_new), lynx_nFem_new(nLynx_new), lynx_rdMortTerr_new(nLynx_new);
      StringVector lynx_breed_new(nLynx_new), lynx_color_new(nLynx_new), lynx_pop_new(nLynx_new), lynx_sex_new(nLynx_new), lynx_status_new(nLynx_new);
      // add residents bit
      for(int i = 0; i<nRes; i++){
        lynx_xcor_new(i) = residents_xcor(i), lynx_ycor_new(i) = residents_ycor(i), lynx_who_new(i) = residents_who(i), lynx_heading_new(i) = residents_heading(i), lynx_prevX_new(i) = residents_prevX(i), lynx_prevY_new(i) = residents_prevY(i), lynx_breed_new(i) = residents_breed(i), lynx_color_new(i) = residents_color(i), lynx_pop_new(i) = residents_pop(i), lynx_sex_new(i) = residents_sex(i), lynx_age_new(i) = residents_age(i), lynx_status_new(i) = residents_status(i), lynx_lastDispX_new(i) = residents_lastDispX(i), lynx_lastDispY_new(i) = residents_lastDispY(i), lynx_nMat_new(i) = residents_nMat(i), lynx_maleID_new(i) = residents_maleID(i), lynx_nFem_new(i) = residents_nFem(i), lynx_rdMortTerr_new(i) = residents_rdMortTerr(i), lynx_steps_new(i) = residents_steps(i);
      }
      
      // List L_return = List::create(Named("int_1") = int_1,
      //                              _["res_index"] = res_index,
      //                              _["nRes"] = nRes,
      //                              _["nDisp_new"] = nDisp_new,
      //                              _["nLynx_new"] = nLynx_new,
      //                              _["dispersers_who_new"] = dispersers_who_new,
      //                              _["dispersers_who"] = dispersers_who,
      //                              _["residents_who"] = residents_who,
      //                              _["ChosenCells_who"] = ChosenCells_who,
      //                              _["lynx_who_new"] = lynx_who_new,
      //                              _["nDispLeft"] = nDispLeft,
      //                              _["step_count"] = step_count);
      // return L_return;
      
      // add res new if any{
      if(nDisp_new>0){
        IntegerVector disp_new_index(nDisp_new);
        for(int i = 0; i<nDisp_new; i++){
          // disp_new_index(i) = i + nRes;
          //}
          // if(nDisp_new < nDisp){
          //   List L_return = List::create(Named("int_100") = int_100,//corresLynx_inter,
          //                                _["step"] = step_count,
          //                                _["nRes"] = nRes,
          //                                _["nDisp"] = nRes,
          //                                _["res_index"] = res_index,
          //                                _["nDisp_new"] = nDisp_new,
          //                                _["deadLynxColl_who"] = deadLynxColl["who"],
          //                                                                    _["deathRoad"] = deathRoad,
          //                                                                    _["disp_new_index"] = disp_new_index,
          //                                                                    _["dispersers_who"] = dispersers_who,
          //                                                                    _["dispersers_steps"] = dispersers_steps,
          //                                                                    _["dispersers_who_new"] = dispersers_who_new,
          //                                                                    _["lynx_who_new"] = lynx_who_new,
          //                                                                    _["lynx_who"] = lynx_who//,
          //   );
          //   return L_return;
          // }
          //lynx_xcor_new[disp_new_index] = dispersers_xcor_new, lynx_steps_new[disp_new_index] = dispersers_steps_new, lynx_ycor_new[disp_new_index] = dispersers_ycor_new, lynx_heading_new[disp_new_index] = dispersers_heading_new, lynx_prevX_new[disp_new_index] = dispersers_prevX_new, lynx_prevY_new[disp_new_index] = dispersers_prevY_new, lynx_breed_new[disp_new_index] = dispersers_breed_new, lynx_color_new[disp_new_index] = dispersers_color_new, lynx_pop_new[disp_new_index] = dispersers_pop_new, lynx_sex_new[disp_new_index] = dispersers_sex_new, lynx_age_new[disp_new_index] = dispersers_age_new, lynx_status_new[disp_new_index] = dispersers_status_new, lynx_maleID_new[disp_new_index] = dispersers_maleID_new, lynx_nFem_new[disp_new_index] = dispersers_nFem_new, lynx_rdMortTerr_new[disp_new_index] = dispersers_rdMortTerr;
          lynx_xcor_new(i + nRes) = dispersers_xcor_new(i), lynx_steps_new(i + nRes) =  dispersers_steps_new(i), lynx_ycor_new(i + nRes) = dispersers_ycor_new(i), lynx_who_new(i + nRes) = dispersers_who_new(i), lynx_heading_new(i + nRes) = dispersers_heading_new(i), lynx_prevX_new(i + nRes) = dispersers_prevX_new(i), lynx_prevY_new(i + nRes) = dispersers_prevY_new(i), lynx_breed_new(i + nRes) = dispersers_breed_new(i), lynx_color_new(i + nRes) = dispersers_color_new(i), lynx_pop_new(i + nRes) = dispersers_pop_new(i), lynx_sex_new(i + nRes) = dispersers_sex_new(i), lynx_age_new(i + nRes) = dispersers_age_new(i), lynx_status_new(i + nRes) = dispersers_status_new(i), lynx_lastDispX_new(i + nRes) = dispersers_lastDispX_new(i), lynx_lastDispY_new(i + nRes) = dispersers_lastDispY_new(i), lynx_nMat_new(i + nRes) = dispersers_nMat_new(i), lynx_maleID_new(i + nRes) = dispersers_maleID_new(i), lynx_nFem_new(i + nRes) = dispersers_nFem_new(i), lynx_rdMortTerr_new(i + nRes) = dispersers_rdMortTerr_new(i);
        }
      }
      // List L_return = List::create(Named("lynx_who_new") = lynx_who_new,
      //                              _["disperser_who_new"] = dispersers_who_new,
      //                              _["lynx_steps_new"] = lynx_steps_new,
      //                              _["disperser_steps_new"] = dispersers_steps_new,
      //                              _["nDisp_new"] = nDisp_new,
      //                              _["deadDisp"] = deadDisp,
      //                              _["deadLynxColl"] = deadLynxColl);
      // return L_return;
      lynx_xcor = clone(lynx_xcor_new), lynx_steps = clone(lynx_steps_new);
      lynx_ycor= clone(lynx_ycor_new), lynx_who= clone(lynx_who_new);
      lynx_heading= clone(lynx_heading_new), lynx_prevX= clone(lynx_prevX_new);
      lynx_prevY= clone(lynx_prevY_new), lynx_breed= clone(lynx_breed_new);
      lynx_color= clone(lynx_color_new), lynx_pop= clone(lynx_pop_new), lynx_sex= clone(lynx_sex_new);
      lynx_age= clone(lynx_age_new), lynx_status= clone(lynx_status_new);
      lynx_lastDispX= clone(lynx_lastDispX_new), lynx_lastDispY= clone(lynx_lastDispY_new); 
      lynx_nMat= clone(lynx_nMat_new), lynx_maleID= clone(lynx_maleID_new); 
      lynx_nFem= clone(lynx_nFem_new), lynx_rdMortTerr= clone(lynx_rdMortTerr_new);
      //}
      
      Rcout << "Rcout 13 lynx updated" << std::endl << step_count << std::endl;
      
      //lynx_xcor = clone(lynx_xcor_new), lynx_steps = clone(lynx_steps_new), 
      //lynx_ycor= clone(lynx_ycor_new), lynx_who= clone(lynx_who_new), 
      //lynx_heading= clone(lynx_heading_new), lynx_prevX= clone(lynx_prevX_new), 
      //lynx_prevY= clone(lynx_prevY_new), lynx_breed= clone(lynx_breed_new), 
      //lynx_color= clone(lynx_color_new), lynx_pop= clone(lynx_pop_new), lynx_sex= clone(lynx_sex_new),
      //lynx_age= clone(lynx_age_new), lynx_status= clone(lynx_status_new), 
      //lynx_lastDispX= clone(lynx_lastDispX_new), lynx_lastDispY= clone(lynx_lastDispY_new), 
      //lynx_nMat= clone(lynx_nMat_new), lynx_maleID= clone(lynx_maleID_new), 
      //lynx_nFem= clone(lynx_nFem_new), lynx_rdMortTerr= clone(lynx_rdMortTerr_new);
      // // List L_return = List::create(Named("lynx_who") = lynx_who,
      // //                              _["lynx_xcor"] = lynx_xcor,
      // //                              _["nDisp_new"] = nDisp_new,
      // //                              _["dispersers_who_new"] = dispersers_who_new,
      // //                              _["deadDisp"] = deadDisp,
      // //                              _["deadLynxColl"] = deadLynxColl
      // //                              );
      // // return L_return;
      
      
      //       
      //       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //       // SEARCHTERRITORY BIT
      
      //compute number females in disp
      int nfemDisp = 0;
      for(int l = 0; l<nDisp_new; l++) {      /* loop over all rows  */
        if(dispersers_sex_new(l) == "F"){
          nfemDisp++;
        }
      }
      //       // if some females do find lines and create subset matrix
      if(nfemDisp > 1){
        IntegerVector index_femDisp(nfemDisp);
        for(int l = 0, j = 0; l<nDisp_new; l++){
          if(dispersers_sex_new(l) == "F"){
            index_femDisp[j] = l;
            j++;
          }
        }
        IntegerVector index_femDispRand = sample(index_femDisp, nfemDisp);
        IntegerVector DispFem_xcor = dispersers_xcor_new[index_femDispRand], DispFem_ycor = dispersers_ycor_new[index_femDispRand], DispFem_who = dispersers_who_new[index_femDispRand], DispFem_heading = dispersers_heading_new[index_femDispRand], DispFem_prevX = dispersers_prevX_new[index_femDispRand], DispFem_prevY = dispersers_prevY_new[index_femDispRand], DispFem_age = dispersers_age_new[index_femDispRand], DispFem_steps = dispersers_steps_new[index_femDispRand], DispFem_lastDispX = dispersers_lastDispX_new[index_femDispRand], DispFem_lastDispY = dispersers_lastDispY_new[index_femDispRand], DispFem_nMat = dispersers_nMat_new[index_femDispRand], DispFem_maleID = dispersers_maleID_new[index_femDispRand], DispFem_nFem = dispersers_nFem_new[index_femDispRand], DispFem_rdMortTerr = dispersers_rdMortTerr_new[index_femDispRand];
        StringVector DispFem_breed = dispersers_breed_new[index_femDispRand], DispFem_color = dispersers_color_new[index_femDispRand], DispFem_pop = dispersers_pop_new[index_femDispRand], DispFem_sex = dispersers_sex_new[index_femDispRand], DispFem_status = dispersers_status_new[index_femDispRand];
        //stop("Got there");
        // List L_return = List::create(Named("DispFem_breed") = DispFem_breed,
        //                              _["nfemDisp"] = nfemDisp,
        //                              _["step"] = step_count
        //                                //_["terrSize"] = terrSize
        // );
        // return L_return;
        //
        // loop for female territory search
        IntegerVector vec_out(nfemDisp);
        for(int f = 0; f<nfemDisp; f++){ //for(searchingFemID in dispFemID) {
          //for(int f = 0; f<2; f++){
          if(
            ((HabitatMap(DispFem_lastDispY(f), DispFem_lastDispX(f)) == int_4) |
              (HabitatMap(DispFem_lastDispY(f), DispFem_lastDispX(f)) == int_3) |
              (HabitatMap(DispFem_lastDispY(f), DispFem_lastDispX(f)) == int_2)) &
              !(R_IsNA(TerrMap(DispFem_lastDispY(f), DispFem_lastDispX(f))))
              //int_0 != int_1
          ){
            int terrSize = 0; // init terrSize
            if((HabitatMap(1, 1) == 4) & R_IsNA(TerrMap(1, 1))){
              for(int l = 0; l<TerrMap.nrow(); l++){
                for(int c = 0; c<TerrMap.ncol(); c++){
                  if(R_IsNA(TerrMap(l,c)) == false){
                    availCellsUpdatedRas(l,c) = int_0;
                    //stop("got in loop");
                  }
                }
              }
              if(popDist(DispFem_lastDispY(f), DispFem_lastDispX(f)) == int_1){
                terrSize = int(coreTerrSizeFAlps);
              }
              if(popDist(DispFem_lastDispY(f), DispFem_lastDispX(f)) == int_2){
                terrSize = int(coreTerrSizeFJura);
              }
              if(popDist(DispFem_lastDispY(f), DispFem_lastDispX(f)) == int_3){
                terrSize = int(coreTerrSizeFVosgesPalatinate);
              }
              if(popDist(DispFem_lastDispY(f), DispFem_lastDispX(f)) == int_4){
                terrSize = int(coreTerrSizeFBlackForest);
              }
            }
            //beginning of spread function /////////////////////////////////////////////////////////////////////////
            
            // some inits
            IntegerMatrix Spredprob = clone(availCellsUpdatedRas);
            IntegerMatrix landscape = clone(availCellsUpdatedRas);
            IntegerVector loci_ind(1); loci_ind(0) = HabitatMap.ncol() * (HabitatMap.nrow() - DispFem_lastDispY(f)) + DispFem_lastDispX(f);
            IntegerVector loci_y(1); loci_y(0) = DispFem_lastDispY(f);
            IntegerVector loci_x(1); loci_x(0) = DispFem_lastDispX(f);
            int MaxSize = terrSize;
            bool ReturnIndices = true;
            bool spreadStateExists = false;
            bool spreadProbLaterExists = true;
            int initialLoci = int(loci_ind(0));
            int sequenceInitialLoci = int_1;
            int ncells = HabitatMap.ncol() * HabitatMap.nrow();
            bool allowOverlapOrReturnDistances = allowOverlap | returnDistances;
            bool useMatrixVersionSpreads = allowOverlapOrReturnDistances | spreadStateExists;
            IntegerVector spreadsDT_spreads(ncells);
            int n = int_1;
            spreadsDT_spreads(loci_ind(0)) = n;
            IntegerVector spreadIndices(100); spreadIndices(0) = 1;
            IntegerVector prevSpreadIndicesActiveLen(1); prevSpreadIndicesActiveLen(0) = 1;
            int prevSpreadIndicesFullLen = spreadIndices.size();
            int size = int_1;
            bool noMaxSize = false;
            ////// now get within while loop
            int iterations = int_5;
            while(n <iterations){
              List potentials = UniqFreeAdjCellsRandOrd(loci_x,
                                                        loci_y,
                                                        availCellsUpdatedRas);
              IntegerVector potentials_AdjX = potentials["AdjX"];
              IntegerVector potentials_AdjY = potentials["AdjY"];
              IntegerVector potentials_CellInd = potentials["CellInd"];
              // potential cells kept
              int nPot = potentials_AdjX.size();
              IntegerVector kept_potentials(0);
              for(int j = 0; j<nPot; j++){
                if(availCellsUpdatedRas(potentials_AdjY(j), potentials_AdjX(j)) == int_1){
                  kept_potentials.push_back(j);
                }//else{
                //kept_potentials(j) = false;
                //}
              }
              IntegerVector potentials_AdjXKept(0);
              IntegerVector potentials_AdjYKept(0);
              IntegerVector potentials_CellIndKept(0);
              if(kept_potentials.size()>int_0){
                for(int i = 0; i<kept_potentials.size();i++){
                  potentials_AdjXKept(i) = potentials_AdjX[kept_potentials(i)];
                  potentials_AdjYKept(i) = potentials_AdjY[kept_potentials(i)];
                  potentials_CellIndKept(i) = potentials_CellInd[kept_potentials(i)];
                }
              }
              int nKeptPot = potentials_CellIndKept.size();// light not be useful for us
              
              if(nKeptPot>0){
                IntegerVector events = clone(potentials_CellIndKept);
                
                if(noMaxSize == false){
                  IntegerVector spreadsDT_spreads_pot = spreadsDT_spreads[potentials_CellIndKept];
                  int len = std::accumulate(spreadsDT_spreads_pot.begin(), spreadsDT_spreads_pot.end(), 0u); // on our case, length is one sone only a sum of that, no need to tabulate  
                  
                  if(((len + size) > MaxSize) & (size < MaxSize)){
                    Rcout << "Rcout spread" << std::endl << step_count << std::endl;
                  }
                  
                }
                
              }
              n++;
            }
            
            ////////////////////////////////////////////////////////////////////////////////////////////////////////
          }// end of habitat 4 spread
        } // for loop over dipersing females
      }// end of dispersin female bit if(n_fem_Disp > 1){
      
      
      /////////////////////////////  
      
    }// of if dispersers left
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // reupdate res and disp tables for next stp
    
    // sizes before table updates
    int nLynx_old = nLynx;
    int nResOld = nRes; // resize does not work in rcpp
    int nDispOld = nDisp; // resize does not work in rcpp
    nDispLeftOld = nDispLeft; // resize does not work in rcpp
    
    // need to kinda reset dispNow as was not carried through above loop and did undergo reduction with death
    nLynx = lynx_status.size();
    if(lynx_dispNow.size() > nLynx){
      lynx_dispNow.erase(nLynx, lynx_dispNow.size());
    }
    
    // if(lynx_lastDispX.size() != lynx_status.size()){
    //   List L_return = List::create(Named("Lynx_status") = lynx_status,
    //                                _["lynx_dispNow"] = lynx_dispNow,
    //                                _["Lynx_lastDispX"] = lynx_lastDispX,
    //                                _["step_count++"] = step_count);
    //   return L_return;
    // }
    
    // now update lynx_dispNow
    for(int l = 0; l<nLynx; l++){
      if((lynx_steps(l) > step) & (lynx_status(l) == "disp")){
        lynx_dispNow(l) = "yes";
      }else{
        lynx_dispNow(l) = "no";
      }
    }
    
    // create disperser and non disperser var /////////////////////////////////////////////////////////////////////
    
    //define some updated dim
    nDisp = N_Eq_Str(lynx_status, "disp");
    nDispLeft = N_Eq_Str(lynx_dispNow, "yes");
    nRes = N_Eq_Str(lynx_dispNow, "no");
    
    Rcout << "Rcout 14 dispnow step" << std::endl << step_count << std::endl;
    // Rcout << "Rcout 14.1 dispnow nDispLeft" << std::endl << nDispLeft << std::endl;
    // Rcout << "Rcout 14.1 dispnow nDispLeftOld" << std::endl << nDispLeftOld << std::endl;
    // Rcout << "Rcout 14.1 dispnow dispersers_who" << std::endl << dispersers_who << std::endl;
    
    // reshape dispersers vector to right size
    if(nDispLeftOld > nDispLeft){
      dispersers_xcor.erase(nDispLeft, nDispLeftOld), dispersers_ycor.erase(nDispLeft, nDispLeftOld), dispersers_who.erase(nDispLeft, nDispLeftOld), dispersers_heading.erase(nDispLeft, nDispLeftOld), dispersers_prevX.erase(nDispLeft, nDispLeftOld), dispersers_prevY.erase(nDispLeft, nDispLeftOld), dispersers_breed.erase(nDispLeft, nDispLeftOld), dispersers_color.erase(nDispLeft, nDispLeftOld), dispersers_pop.erase(nDispLeft, nDispLeftOld), dispersers_sex.erase(nDispLeft, nDispLeftOld), dispersers_age.erase(nDispLeft, nDispLeftOld), dispersers_status.erase(nDispLeft, nDispLeftOld), dispersers_steps.erase(nDispLeft, nDispLeftOld), dispersers_lastDispX.erase(nDispLeft, nDispLeftOld), dispersers_lastDispY.erase(nDispLeft, nDispLeftOld), dispersers_nMat.erase(nDispLeft, nDispLeftOld), dispersers_maleID.erase(nDispLeft, nDispLeftOld), dispersers_nFem.erase(nDispLeft, nDispLeftOld), dispersers_rdMortTerr.erase(nDispLeft, nDispLeftOld);
      // if(dispersers_xcor.size() == nDispLeft){
      //   List L_return = List::create(Named("Where") = "dispersers_xcor.erase",
      //                                _["nDisp"] = nDisp,
      //                                _["Lynx_status"] = lynx_status,
      //                                _["lynx_dispNow"] = lynx_dispNow,
      //                                _["Lynx_lastDispX"] = lynx_lastDispX,
      //                                _["step_count++"] = step_count);
      //   return L_return;
      // }
    }
    
    Rcout << "Rcout 15 dispersers erase update" << std::endl << step_count << std::endl;
    Rcout << "Rcout 15.1 dispersers_who.size()" << std::endl << dispersers_who.size() << std::endl;
    Rcout << "Rcout 15.2 residents_who.size()" << std::endl << residents_steps.size() << std::endl;
    Rcout << "Rcout 15.3 nRes" << std::endl << nRes << std::endl;
    Rcout << "Rcout 15.4 nResOld" << std::endl << nResOld << std::endl;
    Rcout << "Rcout 15.5 lynx_who.size()" << std::endl << lynx_steps.size() << std::endl;
    
    // make loop to fill residents and dispersers vectors again
    for(int i = 0, i_disp = 0, i_ndisp = 0 ; i<nLynx; i++){
      if(lynx_dispNow(i) == "yes"){
        dispersers_xcor(i_disp) = lynx_xcor(i), dispersers_steps(i_disp) = lynx_steps(i),  dispersers_ycor(i_disp) = lynx_ycor(i), dispersers_who(i_disp) = lynx_who(i), dispersers_heading(i_disp) = lynx_heading(i), dispersers_prevX(i_disp) = lynx_prevX(i), dispersers_prevY(i_disp) = lynx_prevY(i), dispersers_breed(i_disp) = lynx_breed(i), dispersers_color(i_disp) = lynx_color(i), dispersers_pop(i_disp) = lynx_pop(i), dispersers_sex(i_disp) = lynx_sex(i), dispersers_age(i_disp) = lynx_age(i), dispersers_status(i_disp) = lynx_status(i), dispersers_lastDispX(i_disp) = lynx_lastDispX(i), dispersers_lastDispY(i_disp) = lynx_lastDispY(i), dispersers_nMat(i_disp) = lynx_nMat(i), dispersers_maleID(i_disp) = lynx_maleID(i), dispersers_nFem(i_disp) = lynx_nFem(i), dispersers_rdMortTerr(i_disp) = lynx_rdMortTerr(i);
        i_disp++;
        // if(i_disp>nDispLeft){
        //   List L_return = List::create(Named("Where") = "dispersers_xcor(i_disp)",
        //                                _["Lynx_status"] = lynx_status,
        //                                _["lynx_dispNow"] = lynx_dispNow,
        //                                _["Lynx_lastDispX"] = lynx_lastDispX,
        //                                _["step_count++"] = step_count);
        //   return L_return;
        // }
      }
      if(lynx_dispNow(i) == "no"){
        if(i_ndisp<nResOld){
          residents_xcor(i_ndisp) = lynx_xcor(i), residents_steps(i_ndisp) = lynx_steps(i),  residents_ycor(i_ndisp) = lynx_ycor(i), residents_who(i_ndisp) = lynx_who(i), residents_heading(i_ndisp) = lynx_heading(i), residents_prevX(i_ndisp) = lynx_prevX(i), residents_prevY(i_ndisp) = lynx_prevY(i), residents_breed(i_ndisp) = lynx_breed(i), residents_color(i_ndisp) = lynx_color(i), residents_pop(i_ndisp) = lynx_pop(i), residents_sex(i_ndisp) = lynx_sex(i), residents_age(i_ndisp) = lynx_age(i), residents_status(i_ndisp) = lynx_status(i), residents_lastDispX(i_ndisp) = lynx_lastDispX(i), residents_lastDispY(i_ndisp) = lynx_lastDispY(i), residents_nMat(i_ndisp) = lynx_nMat(i), residents_maleID(i_ndisp) = lynx_maleID(i), residents_nFem(i_ndisp) = lynx_nFem(i), residents_rdMortTerr(i_ndisp) = lynx_rdMortTerr(i);
          i_ndisp++;
          Rcout << "Rcout 15.61 i_ndisp" << std::endl << i_ndisp << std::endl;
        }else{
          // if(i_ndisp>nRes){
          //   List L_return = List::create(Named("Where") = "residents_xcor(i_ndisp) 2",
          //                                _["Lynx_status"] = lynx_status,
          //                                _["lynx_dispNow"] = lynx_dispNow,
          //                                _["Lynx_lastDispX"] = lynx_lastDispX,
          //                                _["step_count++"] = step_count);
          //   return L_return;
          // }
          residents_xcor.push_back(lynx_xcor(i)), residents_ycor.push_back(lynx_ycor(i)), residents_who.push_back(lynx_who(i)), residents_heading.push_back(lynx_heading(i)), residents_prevX.push_back(lynx_prevX(i)), residents_prevY.push_back(lynx_prevY(i)), residents_breed.push_back(lynx_breed(i)), residents_color.push_back(lynx_color(i)), residents_pop.push_back(lynx_pop(i)), residents_sex.push_back(lynx_sex(i)), residents_age.push_back(lynx_age(i)), residents_status.push_back(lynx_status(i)), residents_steps.push_back(lynx_steps(i)), residents_lastDispX.push_back(lynx_lastDispX(i)), residents_lastDispY.push_back(lynx_lastDispY(i)), residents_nMat.push_back(lynx_nMat(i)), residents_maleID.push_back(lynx_maleID(i)), residents_nFem.push_back(lynx_nFem(i)), residents_rdMortTerr.push_back(lynx_rdMortTerr(i));
          i_ndisp++;
          Rcout << "Rcout 15.62 i_ndisp" << std::endl << i_ndisp << std::endl;
        }
      }
    }
    
    Rcout << "Rcout 16 res and disp updated" << std::endl << step_count << std::endl;
    
  }// end step loop
  
  List L_return = List::create(Named("Lynx_who") = lynx_who,
                               _["lynx_lastDispX"] = lynx_lastDispX,
                               _["lynx_lastDispY"] = lynx_lastDispY,
                               _["dispersers_who"] = dispersers_who,
                               _["residents_who"] = residents_who,
                               _["step_count++"] = step_count);
  //_["terrSize"] = terrSize
  return L_return;
  
  
} // end of function

