#include <Rcpp.h>
using namespace Rcpp;
#include <math.h>       /* atan2 */
#define PI 3.14159265
#include <algorithm>
#include <iostream>
#include <vector>
using namespace std;

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
IntegerVector WhichAbove(IntegerVector ToCheck, int Crit) { // give index cells valua above crit for IntegerVector 
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
IntegerVector WhichEqual(IntegerVector ToCheck, int Crit) { // // give index cells equal for IntegerVector  
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
IntegerVector towards_simple(int x_cur, int y_cur, IntegerVector x_to, IntegerVector y_to){
  IntegerVector degrees(x_to.size());
  for(int i = 0; i<x_to.size(); i++){
    degrees(i) = atan2(y_to(i) - y_cur, x_to(i) - x_cur) * (180 / PI);
    if(degrees(i)<0){
      degrees(i) = degrees(i) + 360;
    }
  }
  return degrees;
}

// [[Rcpp::export]]
// towards simple,  convert radians to degrees
int towards_simple_unique(int x_cur, int y_cur, int x_to, int y_to){
  int degrees = 0;
  degrees = atan2(y_to - y_cur, x_to - x_cur) * (180 / PI);
  if(degrees<0){
    degrees = degrees + 360;
  }
  return degrees;
}

// [[Rcpp::export]]
// change heading for some individual whom steps have been reset
int changeHeading(int to_change){
  IntegerVector headings = {0, 45, 90, 135, 180, 225, 270, 315};
  IntegerVector deltas(headings.size());
  for(int i = 0; i<deltas.size(); i++){
    deltas(i) = abs(headings(i) - to_change);
  }
  int delta_min = *min_element(deltas.begin(), deltas.end());
  IntegerVector new_val = WhichEqual(deltas, delta_min);
  return headings(new_val(0));
}

// [[Rcpp::export]]
// adjacent cells coordinates
List UniqAdjCells(IntegerVector x_coords, IntegerVector y_coords, NumericMatrix Matrix){
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
        if((new_y>=0) & (new_y<nRowMat) & (new_x>=0) & (new_x<nColMat)){
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
  IntegerVector AdjX_left = AdjX[keptCells];
  IntegerVector AdjY_left = AdjY[keptCells];
  IntegerVector CellInd_left = CellInd[keptCells];
  List L_return = List::create(Named("AdjX") = AdjX_left,
                          _["AdjY"] = AdjY_left,
                          _["CellInd"] = CellInd_left);
  // return L_return;
  return L_return;
}

// [[Rcpp::export]]
// shorten vector
IntegerVector ShortenIntVec(IntegerVector LongVec, int newSize){
int oldSize = LongVec.size();
  if(oldSize>newSize){
    LongVec.erase(newSize, oldSize);
  }
  return LongVec;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

///*** R
//timesTwo(42)
//*/
