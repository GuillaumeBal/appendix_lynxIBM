#include <Rcpp.h>
using namespace Rcpp;
#include <math.h>       /* atan2 */
#define PI 3.14159265
#include <cstdlib>
#include <iostream>
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
// here I simply explore dynamic allocation for when size unknown in advance

// [[Rcpp::export]]
IntegerVector dynaVec(IntegerVector Length_IntVec,
                      int Length_Subset){
  
  std::vector<int> output;
  //for(int i = 0; i<Length_Subset; i++){
    output = {sample(Length_IntVec(), Length_Subset, false)}; //number of samples is argument 2
  //}
  return Rcpp::IntegerVector(output.begin(), output.end());
  
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

///*** R
//timesTwo(42)
//*/
