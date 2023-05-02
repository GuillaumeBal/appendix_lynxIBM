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
List spreadGB(// DataFrame, NumericVector
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
    bool allowOverlap,
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
  int int_100 = 100;
  
  //beginning of spread function /////////////////////////////////////////////////////////////////////////
  
  // some inits
  IntegerMatrix Spredprob = clone(availCellsUpdatedRas);
  IntegerMatrix landscape = clone(availCellsUpdatedRas);
  IntegerVector loci_ind(1); loci_ind(0) = HabitatMap.ncol() * (HabitatMap.nrow() - DispFem_lastDispY) + DispFem_lastDispX;
  IntegerVector loci_y(1); loci_y(0) = DispFem_lastDispY;
  IntegerVector loci_x(1); loci_x(0) = DispFem_lastDispX;
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
  // sanity check var def
  Rcout << "Rcout spread 1 var def prevSpreadIndicesFullLen : " << std::endl << prevSpreadIndicesFullLen << std::endl;
  //// now get within while loop
  int iterations = int_5;
  while((n - 1) < iterations){
    n++;
    // sanity check
    Rcout << "Rcout spread 4 n : " << std::endl << n << std::endl;
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
    // sanity check
    //Rcout << "Rcout spread 2 / kept_potentials" << std::endl <<  kept_potentials << std::endl;
    
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
    
    // sanity check
    //Rcout << "Rcout spread 3 / nKeptPot" << std::endl <<  nKeptPot << std::endl;
    
    if(nKeptPot>0){
      IntegerVector events = clone(potentials_CellIndKept);
      
      if(noMaxSize == false){
        IntegerVector spreadsDT_spreads_pot = spreadsDT_spreads[potentials_CellIndKept];
        int len = std::accumulate(spreadsDT_spreads_pot.begin(), spreadsDT_spreads_pot.end(), 0u); // on our case, length is one sone only a sum of that, no need to tabulate
        
        if(((len + size) > MaxSize) & (size < MaxSize)){
          // sanity check
          //Rcout << "Rcout spread 3 / len" << std::endl << len << std::endl;
          int stuff = 0;
        }
        
      }
      
    }
    
  }// while(n< iterations)
  
  List L_return = List::create(Named("where") = "very end",
                               _["prevSpreadIndicesFullLen"] = prevSpreadIndicesFullLen,
                               _["size"] = prevSpreadIndicesFullLen);
  return L_return;
  
} // end of function