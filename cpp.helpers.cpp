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
IntegerVector IntVecSubIndex(IntegerVector ToSub, IntegerVector PosToKeep) { // check number cells vector matching character string 
 IntegerVector kept(PosToKeep.size());
   for(int i = 0; i<PosToKeep.size(); i++){
     kept(i) = ToSub(PosToKeep(i));
   }
  return kept;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

///*** R
//timesTwo(42)
//*/
