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
  IntegerVector which_vec_above(n_above);
  int p = 0;
  for(int i = 0; i<ToCheck.size(); i++){
    if(ToCheck(i)>Crit){
      which_vec_above(p) = i;
      p++;
    }
  }
  return which_vec_above;
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
  IntegerVector which_vec_equal(n_equal);
  int p = 0;
  for(int i = 0; i<ToCheck.size(); i++){
    if(ToCheck(i)==Crit){
      which_vec_equal(p) = i;
      p++;
    }
  }
  return which_vec_equal;
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
  IntegerVector which_vec_in_set(n_in);
  int p = 0;
  for(int i = 0; i<ToCheck.size(); i++){
    for(int j = 0; j<subset.size(); j++){
      if(ToCheck(i)==subset(j)){
        which_vec_in_set(p) = i;
        p++;
      }
    }
  }
  return which_vec_in_set;
}

// [[Rcpp::export]]
// check cells within subset for integer 
bool IntInSet(int ToCheck, IntegerVector pool) { 
  bool is_in_set = false;
  for(int p = 0; p<pool.size(); p++){
    if(pool(p) == ToCheck){
      is_in_set = true;
      p = pool.size() - 1;
    }
  }
  return is_in_set;
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
// cell number to coords on map in cpp
List CellNumtoRowCol(IntegerVector cellNum, IntegerMatrix Matrix){
  IntegerVector rowNum(cellNum.size());
  IntegerVector colNum(cellNum.size());
  int nRowMat = Matrix.nrow();
  int nColMat = Matrix.ncol();
  for(int i = 0; i<cellNum.size(); i++){
    rowNum(i) = nRowMat - trunc(cellNum(i) / nColMat) - 1;
    //trunc(cellNum(i) / nColMat);
    colNum(i) = cellNum(i) - nColMat * (nRowMat - rowNum(i) - 1) - 1 ; 
    //cellNum(i) - (rowNum(i) * nColMat + 1);
  }
  List coordsRC = List::create(Named("cellNum") = cellNum,
                               _["x_coords"] = colNum,
                               _["y_coords"] = rowNum);
  return coordsRC;
}

// [[Rcpp::export]]
// x y to cellNum
IntegerVector RowColtoCellNum(IntegerVector x_coords, IntegerVector y_coords, IntegerMatrix Matrix){
  IntegerVector cellNum(y_coords.size());
  int nRowMat = Matrix.nrow();
  int nColMat = Matrix.ncol();
  for(int i = 0; i<cellNum.size(); i++){
    cellNum(i) = nColMat * (nRowMat - (y_coords(i) + 1)) + (x_coords(i) + 1);
    //HabitatMap.ncol() * (HabitatMap.nrow() - (DispFem_lastDispY + 1)) + (DispFem_lastDispX + 1);
  }
  return cellNum;
}

// [[Rcpp::export]]
// adjacent cells coordinates
List UniqFreeAdjCellsRandOrd(IntegerVector cellNum, IntegerMatrix Matrix){
  IntegerVector deltaX = {-1, 0, 1};
  IntegerVector deltaY = {-1, 0, 1};
  int nColMat = Matrix.ncol();
  int nRowMat = Matrix.nrow();
  int int_1 = 1;
  IntegerVector AdjX(0);
  IntegerVector AdjY(0);
  IntegerVector CellInd(0);
  int new_y;
  int new_x;
  int new_index;
  // change cell num into XY
  List all_coords = CellNumtoRowCol(cellNum = cellNum, Matrix = Matrix);
  IntegerVector x_coords = all_coords["x_coords"];
  IntegerVector y_coords = all_coords["y_coords"];
  //Rcout << "Rcout 1 inits" << std::endl << CellInd << std::endl;
  for(int z = 0; z<x_coords.size(); z++){
    //Rcout << "Rcout z" << std::endl << z << std::endl;
    for(int l = 0; l<deltaY.size(); l++){
      for(int c = 0; c<deltaX.size(); c++){
        new_y = y_coords(z) + deltaY(l);
        new_x = x_coords(z) + deltaX(c);
        //HabitatMap.ncol() * (HabitatMap.nrow() - DispFem_lastDispY) + DispFem_lastDispX;
        //ncol(my.mat) * (nrow(my.mat) - my.coords[1]) + my.coords[2]
        new_index = nColMat * (nRowMat -  (new_y + 1)) + (new_x + 1); // plus one because starts 0 in cpp
        bool new_index_exist = IntInSet(new_index, cellNum); // check does not already exist
        if((new_y>=0) & (new_y<nRowMat) & (new_x>=0) & (new_x<nColMat) & 
           new_index_exist == false &
           //availCellsUpdatedRas.nrow() - 1 - potentials_AdjY(j)
           Matrix(nRowMat - 1 - new_y, new_x) == int_1){ // add to keep only unoccupied cells
          AdjX.push_back(new_x);
          AdjY.push_back(new_y);
          CellInd.push_back(new_index);
        }
      }
    }
  }
  //Rcout << "Rcout cells with duplicates" << std::endl << CellInd << std::endl;
  // find unique index and then pick just one of each
  IntegerVector UniqCell = unique(CellInd);
  //Rcout << "Rcout UniqCell" << std::endl << UniqCell << std::endl;
  int nUniqCell = UniqCell.size();
  IntegerVector keptCells(nUniqCell);
  for(int i = 0; i<nUniqCell; i++){
    for(int j = 0; j<CellInd.size(); j++){
      if(CellInd(j) == UniqCell(i)){
        keptCells(i) = j; // this is a position
        //break;
      }
    }
  }
  //Rcout << "Rcout keptCells" << std::endl << keptCells << std::endl;
  IntegerVector randorder = sample(keptCells.size(), keptCells.size(), false) - 1;
  IntegerVector keptCellsRandOrder = keptCells[randorder];
  IntegerVector AdjX_left = AdjX[keptCellsRandOrder];
  IntegerVector AdjY_left = AdjY[keptCellsRandOrder];
  IntegerVector CellInd_left = CellInd[keptCellsRandOrder];
  //Rcout << "Rcout CellInd_left" << std::endl << CellInd_left << std::endl;
  List L_return = List::create(Named("CellInd") = CellInd_left,
                               _["AdjX"] = AdjX_left,
                               _["AdjY"] = AdjY_left);
  // return L_return;
  return L_return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
// the spread function, fire like
List spreadGB(// DataFrame, NumericVector
    IntegerMatrix HabitatMap, // different from availCellUpdateRas ?
    IntegerMatrix availCellsUpdatedRas,
    int DispFem_lastDispY,
    int DispFem_lastDispX,
    int terrSize
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
  int int_10 = 10;
  int int_20 = 20;
  int int_100 = 100;
  int int_1e3 = 1000;
  int int_1e4 = 10000;
  int int_1e5 = 100000;
  int int_1e6 = 1000000;
  
  //beginning of spread function /////////////////////////////////////////////////////////////////////////
  
  // some inits
  IntegerVector loci_ind(1); 
  //nColMat * (nRowMat -  (new_y + 1)) + (new_x + 1)
  loci_ind(0) = HabitatMap.ncol() * (HabitatMap.nrow() - (DispFem_lastDispY + 1)) + (DispFem_lastDispX + 1);
  IntegerVector loci_y(1); loci_y(0) = DispFem_lastDispY;
  IntegerVector loci_x(1); loci_x(0) = DispFem_lastDispX;
  int MaxSize = terrSize;
  int initialLoci = int(loci_ind(0));
  int sequenceInitialLoci = int_1;
  int nCells = HabitatMap.ncol() * HabitatMap.nrow();
  int n = int_1;
  
  //// now get within while loop
  int iterations =  int_20;//int_1e6;
  IntegerVector events; // define high because used everywhere 
  
  ////////////////////////////////////////////////////////////////////////////////////
  // while loop
  
  while((loci_ind.size() < MaxSize) & ((n - 1) < iterations)){ 
    
    ///////////////////////////////////////////////
    //find potential cells
    
    // sanity check
    // Rcout << "Rcout spread 4 n : " << std::endl << n << std::endl;
    
    
    //loci_x and loci_y, before when was based on XY coord
    //Rcout << "Rcout loci_ind start while: " << std::endl << loci_ind << std::endl;
    List potentials = UniqFreeAdjCellsRandOrd(loci_ind, availCellsUpdatedRas);
    // sanity check
    //return(potentials);
    IntegerVector potentials_AdjX = potentials["AdjX"];
    IntegerVector potentials_AdjY = potentials["AdjY"];
    IntegerVector potentials_CellInd = potentials["CellInd"];
    Rcout << "Rcout potentials_CellInd: " << std::endl << potentials_CellInd << std::endl;
    // potential cells kept
    int nPot = potentials_AdjX.size();
    
    if(nPot > 0){
      // decided to put everything in a single loop
      for(int p = 0; p<nPot; p++){
        if(loci_ind.size() < MaxSize){
          loci_ind.push_back(potentials_CellInd(p));
        }else{
          n = iterations;
        }
      }
      //Rcout << "Rcout loci_ind updated: " << std::endl << loci_ind << std::endl;
    }else{ //if not_a_closure new potential cells
      n = iterations;
    }
    
    n++;// why here ?
    //Rcout << "Rcout n value update : " << std::endl << n << std::endl;
    
  }// while(n< iterations)
  
  //loci_ind.erase(loci_ind.begin());// # in fact keep first cell that is female location
  
  //Rcout << "Rcout n" << std::endl << n << std::endl;
  List L_return_4 = List::create(Named("initialLoci") = initialLoci,
                                 _["cell_ind"] = loci_ind 
  );
  return L_return_4;

}// end of function
