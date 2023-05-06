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
  //Rcout << "Rcout 1 inits" << std::endl << CellInd << std::endl;
  for(int z = 0; z<x_coords.size(); z++){
    for(int l = 0; l<deltaY.size(); l++){
      for(int c = 0; c<deltaX.size(); c++){
        new_y = y_coords(z) + deltaY(l);
        new_x = x_coords(z) + deltaX(c);
        //HabitatMap.ncol() * (HabitatMap.nrow() - DispFem_lastDispY) + DispFem_lastDispX;
        //ncol(my.mat) * (nrow(my.mat) - my.coords[1]) + my.coords[2]
        new_index = nColMat * (nRowMat -  (new_y + 1)) + (new_x + 1); // plus one because starts 0 in cpp
        if((new_y>=0) & (new_y<nRowMat) & (new_x>=0) & (new_x<nColMat) &
           ((deltaY(l) != 0) | (deltaX(c) != 0))){
          if(Matrix(new_y, new_x) == 0 ){ // add to keep only unoccupied cell  s
            AdjX.push_back(new_x);
            AdjY.push_back(new_y);
            CellInd.push_back(new_index);
            //stop("Went into loop");
            //if(new_index >(nColMat * nRowMat)){
            //  Rcout << "Index too big" << std::endl << new_index << std::endl;
            //}
          }
        }
      }
    }
  }
  //Rcout << "Rcout 2 cell with duplicates" << std::endl << CellInd << std::endl;
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
  //Rcout << "Rcout 3 cells kept" << std::endl << keptCells << std::endl;
  IntegerVector randorder = sample(AdjX.size(), AdjX.size(), false) - 1;
  IntegerVector keptCellsRandOrder = keptCells[randorder];
  IntegerVector AdjX_left = AdjX[keptCellsRandOrder];
  IntegerVector AdjY_left = AdjY[keptCellsRandOrder];
  IntegerVector CellInd_left = CellInd[keptCellsRandOrder];
  //Rcout << "Rcout 4 before outputs" << std::endl << keptCells << std::endl;
  List L_cells = List::create(Named("AdjX") = AdjX_left,
                              _["AdjY"] = AdjY_left,
                              _["CellInd"] = CellInd_left);
  // return L_return;
  return L_cells;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
// the spread function, fire like
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
  int int_1e3 = 1000;
  int int_1e4 = 10000;
  int int_1e5 = 100000;
  int int_1e6 = 1000000;
  
  //beginning of spread function /////////////////////////////////////////////////////////////////////////
  
  // some inits
  IntegerMatrix Spredprob = clone(availCellsUpdatedRas);
  IntegerMatrix landscape = clone(availCellsUpdatedRas);
  IntegerVector loci_ind(1); 
  //nColMat * (nRowMat -  (new_y + 1)) + (new_x + 1)
  loci_ind(0) = HabitatMap.ncol() * (HabitatMap.nrow() - (DispFem_lastDispY + 1)) + (DispFem_lastDispX + 1);
  IntegerVector loci_y(1); loci_y(0) = DispFem_lastDispY;
  IntegerVector loci_x(1); loci_x(0) = DispFem_lastDispX;
  int MaxSize = terrSize;
  bool ReturnIndices = true;
  bool spreadStateExists = false;
  bool spreadProbLaterExists = true;
  int initialLoci = int(loci_ind(0));
  int sequenceInitialLoci = int_1;
  int nCells = HabitatMap.ncol() * HabitatMap.nrow();
  bool allowOverlapOrReturnDistances = allowOverlap | returnDistances;
  bool useMatrixVersionSpreads = allowOverlapOrReturnDistances | spreadStateExists;
  IntegerVector spreadsDT_spreads(nCells);
  int n = int_1;
  spreadsDT_spreads(initialLoci - 1) = n; // as cpp starts at 0
  IntegerVector spreadIndices(100); spreadIndices(0) = 1;
  IntegerVector prevSpreadIndicesActiveLen(1); prevSpreadIndicesActiveLen(0) = 1;
  int prevSpreadIndicesFullLen = spreadIndices.size();
  int size = int_1;
  bool noMaxSize = false;
  // sanity check var def
  //Rcout << "Rcout 1.1 loci : " << std::endl << loci_ind << std::endl;
  //Rcout << "Rcout 1.2 nCells : " << std::endl << nCells << std::endl;
  //Rcout << "Rcout 1.3 prevSpreadIndicesFullLen : " << std::endl << prevSpreadIndicesFullLen << std::endl;
  //Rcout << "Rcout 1.4 spreadsDT_spreads : " << std::endl << spreadsDT_spreads << std::endl;
  //Rcout << "Rcout 1.5 spreadIndices : " << std::endl << spreadIndices << std::endl;
  //stop("beginning check");
  // List L_return_1 = List::create(Named("where") = "base def parameters",
  //                              _["loci_ind"] = loci_ind,
  //                              _["prevSpreadIndicesFullLen"] = prevSpreadIndicesFullLen,
  //                              _["spreadsDT_spreads"] = spreadsDT_spreads,
  //                              _["spreadIndices"] = spreadIndices);
  // return L_return_1;
  
  //// now get within while loop
  int iterations = int_1e6;
  IntegerVector events; // define high because used everywhere 
  
  ////////////////////////////////////////////////////////////////////////////////////
  // while loop
  
  while((loci_ind.size() <= int_1) & ((n - 1) < iterations)){
    
    ///////////////////////////////////////////////
    //find potential cells
    
    // sanity check
    // Rcout << "Rcout spread 4 n : " << std::endl << n << std::endl;
    List potentials = UniqFreeAdjCellsRandOrd(loci_x,
                                              loci_y,
                                              availCellsUpdatedRas);
    // // sanity check
    // return(potentials);
    
    IntegerVector potentials_AdjX = potentials["AdjX"];
    IntegerVector potentials_AdjY = potentials["AdjY"];
    IntegerVector potentials_CellInd = potentials["CellInd"];
    // potential cells kept
    int nPot = potentials_AdjX.size();
    IntegerVector kept_potentials(0);
    if(nPot>0){
      for(int j = 0; j<nPot; j++){
        // if(((availCellsUpdatedRas.nrow() - 1) < potentials_AdjY(j)) | ((availCellsUpdatedRas.ncol() - 1) < potentials_AdjX(j))){
        //   stop("you fucked up, too big");
        // }
        // if((0 > potentials_AdjY(j)) | (0 > potentials_AdjX(j))){
        //   stop("you fucked up, too small");
        // }
        if(availCellsUpdatedRas(potentials_AdjY(j), potentials_AdjX(j)) == int_0){
          kept_potentials.push_back(j);
        }//else{
        //kept_potentials(j) = false;
        //}
      }
    }
    int nKeptPot = kept_potentials.size();
    // sanity check
    //Rcout << "Rcout spread 2 / kept_potentials" << std::endl <<  nKeptPot << std::endl;
    // List L_return_2 = List::create(Named("where") = "potentials kept",
    //                                _["loci_ind"] = loci_ind,
    //                                _["potentials_CellInd"] = potentials_CellInd,
    //                                _["kept_potentials"] = kept_potentials,
    //                                _["nKeptPot"] = nKeptPot);
    // return L_return_2;
    
    IntegerVector potentials_AdjXKept(nKeptPot);
    IntegerVector potentials_AdjYKept(nKeptPot);
    IntegerVector potentials_CellIndKept(nKeptPot);
    
    n++;// why here ?
    
    if(nKeptPot>int_0){
      for(int i = 0; i<nKeptPot;i++){
        potentials_AdjXKept(i) = potentials_AdjX[kept_potentials(i)];
        potentials_AdjYKept(i) = potentials_AdjY[kept_potentials(i)];
        potentials_CellIndKept(i) = potentials_CellInd[kept_potentials(i)];
      }
    }
    
    // sanity check
    //Rcout << "Rcout spread 3 / nKeptPot" << std::endl <<  nKeptPot << std::endl;
    // List L_return_3 = List::create(Named("where") = "potentials kept",
    //                                _["loci_ind"] = loci_ind,
    //                                _["potentials_CellInd"] = potentials_CellInd,
    //                                _["kept_potentials"] = kept_potentials,
    //                                _["potentials_CellIndKept"] = potentials_CellIndKept);
    // return L_return_3;
    
    
    ///////////////////////////////////////
    // if som potential cells
    
    if(nKeptPot>0){
      
      events = clone(potentials_CellIndKept);
      
      if(noMaxSize == false){
        IntegerVector spreadsDT_spreads_pot = spreadsDT_spreads[potentials_CellIndKept];
        int len = std::accumulate(spreadsDT_spreads_pot.begin(), spreadsDT_spreads_pot.end(), int_0); // on our case, length is one sone only a sum of that, no need to tabulate
        
        if(((len + size) > MaxSize) & (size < MaxSize)){
          // sanity check
          //Rcout << "Rcout spread 3 / len" << std::endl << len << std::endl;
          // numer of active cells to remove to no more than size of territory possible
          //int toRm = (size + len) - MaxSize; replace by to keep
          IntegerVector cells_to_kep = sample(spreadsDT_spreads_pot.size(), MaxSize, false) - 1;//number of samples is argument 2, -1 to make it start at 0
          IntegerVector spreadsDT_spreads_pot_resized = spreadsDT_spreads_pot[cells_to_kep];
          spreadsDT_spreads_pot = clone(spreadsDT_spreads_pot_resized);
          events = clone(spreadsDT_spreads_pot);
        }
        //Rcout << "Rcout len" << std::endl << len << std::endl;
        //Rcout << "Rcout len, n" << std::endl << n << std::endl;
        size = std::min((size + len) , MaxSize);
      } // if(noMaxSize == false){
      
      // sanity check
      // if(n > int_1e3){
      // List L_return_4 = List::create(Named("where") = "after noMaxSize loop",
      //                                _["loci_ind"] = loci_ind,
      //                                _["events"] = events,
      //                                _["size"] = size);
      // return L_return_4;
      // }
      
      // what to do depending of length of events
      if(events.size()>0){
        int curEventsLen = events.size();
        IntegerVector addedIndices(curEventsLen);
        for(int i = 0; i<curEventsLen; i++){
          // maybe error and  prevSpreadIndicesActiveLen shoud be integer
          addedIndices(i) = prevSpreadIndicesActiveLen(0) + i + int_1;
        }
        // if more values push back
        
        if((curEventsLen + prevSpreadIndicesActiveLen(0)) > prevSpreadIndicesFullLen){
          for(int j= prevSpreadIndicesFullLen; j<(curEventsLen + prevSpreadIndicesActiveLen(0)); j++){
            spreadIndices.push_back(0); //increase length if necessary, then complete below 
          }
        }
        for(int j = 0; j<addedIndices.size(); j++){
          spreadIndices(addedIndices(j)) = events(j);
        }
      } //if(events.size()>0){
      
    }else{// if(nKeptPot>0){
      IntegerVector events_null;
      events = clone(events_null);
    }
    IntegerVector loci_ind_null;
    loci_ind = clone(events);
    //loci_ind = clone(loci_ind_null);
    //loci_ind = clone(loci_ind_null);
    //loci_y set to null as well ?
    //loci_x
  }// while(n< iterations)
  
  Rcout << "Rcout n" << std::endl << n << std::endl;
  // List L_return_4 = List::create(Named("where") = "after while loop",
  //                                _["prevSpreadIndicesActiveLen"] = prevSpreadIndicesActiveLen,
  //                                _["events"] = events,
  //                                _["size"] = size);
  // return L_return_4;
  
  /////////////////////////////////////////////////////////////////////////
  // after while loop
  IntegerVector spreadsIndices_final_indices(prevSpreadIndicesActiveLen(0)); 
  for(int j = 0; j < prevSpreadIndicesActiveLen(0); j++){
    spreadsIndices_final_indices(j) = j;
  }
  IntegerVector spreadIndices_final = spreadIndices[spreadsIndices_final_indices]; 
  IntegerVector initialLocus_final(spreadIndices_final.size()); 
  StringVector id_final(spreadIndices_final.size());
  LogicalVector active_final(spreadIndices_final.size());
  for(int j = 0; j < prevSpreadIndicesActiveLen(0); j++){
    initialLocus_final(j) = initialLoci;
    id_final(j) = int_1;
    active_final(j) = false;
  }
  
  List L_return = List::create(Named("where") = "very end",
                               _["initialLocus_final"] =  initialLocus_final,
                               _["id_final"] = id_final,
                               _["active_final"] = active_final);
  return L_return;
  
}// end of function
