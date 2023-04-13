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




// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

///*** R
//timesTwo(42)
//*/
