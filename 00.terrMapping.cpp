#include <Rcpp.h>
using namespace Rcpp;
#include <math.h>       /* atan2 */
#define PI 3.14159265
#include <algorithm>


////////////////////////////////////////////////////////////////////////////////

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
// cell number to coords on map in cpp
List CellNumtoRowCol(IntegerVector cellNum, IntegerMatrix Matrix, LogicalVector cpp_as_input, LogicalVector cpp_as_output){
  IntegerVector rowNum_xy(cellNum.size());
  IntegerVector colNum_xy(cellNum.size());
  int nRowMat = Matrix.nrow();
  int nColMat = Matrix.ncol();
  int base_rowNum = -100;
  IntegerVector cellNumRC(cellNum.size()); 
  for(int i = 0; i<cellNum.size(); i++){
    if(cpp_as_input(0) == FALSE){
      cellNumRC(i) = cellNum(i);
    }else{
      //Rcout << "Rcout cpp_as_input cell bef " << std::endl << cellNum(i) << std::endl;
      cellNumRC(i) = cellNum(i) + 1;
      //Rcout << "Rcout cpp_as_input " << std::endl << cpp_as_input(0) << std::endl;
      //Rcout << "Rcout cpp_as_input cell after " << std::endl << cellNum(i) << std::endl;
    }
    base_rowNum =  round(cellNumRC(i) / nColMat) ;
    rowNum_xy(i) = (nRowMat - base_rowNum + 1); 
    colNum_xy(i) = (cellNumRC(i) - (base_rowNum - 1) * nColMat); // same, correction to compute regular, and then substract 1
    if(colNum_xy(i) > nColMat){ // issue with constant rounding down whatever function used
      base_rowNum = base_rowNum + 1;
      rowNum_xy(i) = (nRowMat - base_rowNum + 1); // correction to make computations with regular space and then some substraction at the end for c++ scale
      colNum_xy(i) = (cellNumRC(i) - (base_rowNum - 1) * nColMat); // same, correction to compute regular, and then substract 1
    }
    if(cpp_as_output(0) == TRUE){// is boolean
      rowNum_xy(i) = rowNum_xy(i) - 1 ; // correction to make computations with regular space and then some substraction at the end for c++ scale
      colNum_xy(i) = colNum_xy(i) - 1 ;
      cellNumRC(i) = cellNumRC(i) - 1;
      //Rcout << "Rcout cpp_as_output " << std::endl << cpp_as_output(0) << std::endl;
    }
  }
  List coordsRC = List::create(Named("cellNum") = cellNumRC,
                               //_["base_rowNum"] = base_rowNum,
                               _["x_coords"] = colNum_xy,
                               _["y_coords"] = rowNum_xy);
  return coordsRC;
}

// [[Rcpp::export]]
// adjacent cells coordinates
List UniqFreeAdjCellsRandOrd(IntegerVector cellNum, IntegerMatrix Matrix, LogicalVector cpp_as_input, LogicalVector cpp_as_output){
  IntegerVector deltaX = {-1, 0, 1};
  IntegerVector deltaY = {-1, 0, 1};
  int nColMat = Matrix.ncol();
  int nRowMat = Matrix.nrow();
  int int_1 = 1;
  int int_2 = 2;
  IntegerVector AdjX(0);
  IntegerVector AdjY(0);
  IntegerVector cellIndex(0);
  int new_y;
  int new_x;
  int new_index;
  if(cpp_as_input(0) == 0){
    for(int i = 0; i < cellNum.size(); i++){
      cellNum(i) = cellNum(i) - 1;
    }
  }
  // change cell num into XY
  //List L_1 = List::create(Named("int_1") = int_1);
  //return(L_1);
  //Rcout << "Rcout before cellIndex " << std::endl << int_1 << std::endl;
  List all_coords = CellNumtoRowCol(cellNum = cellIndex, Matrix = Matrix, cpp_as_input = 0, cpp_as_output = 0);
  //List L_2 = List::create(Named("all_coords") = all_coords);
  //return(L_2);
  IntegerVector x_coords = all_coords["x_coords"];
  IntegerVector y_coords = all_coords["y_coords"];
  //Rcout << "Rcout 1 inits" << std::endl << cellIndex << std::endl;
  for(int z = 0; z<x_coords.size(); z++){
    //Rcout << "Rcout z" << std::endl << z << std::endl;
    for(int l = 0; l<deltaY.size(); l++){
      for(int c = 0; c<deltaX.size(); c++){
        new_y = y_coords(z) + deltaY(l);
        new_x = x_coords(z) + deltaX(c);
        //HabitatMap.ncol() * (HabitatMap.nrow() - YCoordinate) +XCoordinate;
        //ncol(my.mat) * (nrow(my.mat) - my.coords[1]) + my.coords[2]
        new_index = nColMat * (nRowMat -  (new_y + 1)) + (new_x + 1); // plus one because starts 0 in cpp
        // if(new_index<0){
        //List L_3 = List::create(Named("new_index") = new_index);
        //return(L_3);
        // }
        bool new_index_exist = IntInSet(new_index, cellIndex); // check does not already exist
        if(
          (new_index>=0) &
            (new_y>=0) & (new_y<nRowMat) & (new_x>=0) & (new_x<nColMat) & 
            (new_index_exist == false)){ // add to keep only unoccupied cells
          if((Matrix(nRowMat - 1 - new_y, new_x) == int_1)){
            AdjX.push_back(new_x);
            AdjY.push_back(new_y);
            cellIndex.push_back(new_index);
          }
        }
      }
    }
  }
  //Rcout << "Rcout cells with duplicates" << std::endl << cellIndex << std::endl;
  // find unique index and then pick just one of each
  IntegerVector UniqCell = unique(cellIndex);
  //Rcout << "Rcout UniqCell" << std::endl << UniqCell << std::endl;
  int nUniqCell = UniqCell.size();
  IntegerVector keptCells(nUniqCell);
  for(int i = 0; i<nUniqCell; i++){
    for(int j = 0; j<cellIndex.size(); j++){
      if(cellIndex(j) == UniqCell(i)){
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
  IntegerVector cellIndex_left = cellIndex[keptCellsRandOrder];
  //Rcout << "Rcout cellIndex_left" << std::endl << cellIndex_left << std::endl;
  if(cpp_as_output(0) == 0){
    for(int i = 0; i < cellIndex_left.size(); i++){
      cellIndex_left(i) = cellIndex_left(i) + 1;
      AdjX_left(i) = AdjX_left(i) + 1;
      AdjY_left(i) = AdjY_left(i) + 1;
    }
  }
  List L_return = List::create(Named("CellNum") = cellIndex_left,
                               _["AdjX"] = AdjX_left,
                               _["AdjY"] = AdjY_left);
  // return L_return;
  return L_return;
}

////////////////////////////////////////////////////////////////////////////////
// spread function

// [[Rcpp::export]]
// the spread function, fire like
List spreadGB(// DataFrame, NumericVector
    IntegerMatrix availCellsMat,
    int YCoordinate,
    int XCoordinate,
    int terrSizeMax,
    LogicalVector cpp_as_input,
    LogicalVector cpp_as_output
){
  
  ///////////////////////////////////////////////////////////////////////////////////////
  // TRICKS AND FUNCTIONS FOR CODE
  
  if(cpp_as_input(0) == 0){
    XCoordinate = XCoordinate - 1;
    YCoordinate = YCoordinate - 1;
  }
  
  IntegerMatrix HabitatMap = clone(availCellsMat);
  
  // some integer def because IntegerVector v(1) = 1 does not work
  int int_0 = 0;
  int int_1 = 1;
  int int_5 = 5;
  int int_10 = 10;
  int int_100 = 100;
  
  //beginning of spread function /////////////////////////////////////////////////////////////////////////
  
  // some inits
  IntegerVector loci_ind(1); 
  //nColMat * (nRowMat -  (new_y + 1)) + (new_x + 1)
  loci_ind(0) = (HabitatMap.nrow() - 
    (YCoordinate + 1)) *  // would be + 0 if not cpp 
    HabitatMap.ncol() + XCoordinate + 1; // plus 1 because Xcoord cpp 
  IntegerVector loci_y(1); loci_y(0) = YCoordinate;
  IntegerVector loci_x(1); loci_x(0) = XCoordinate;
  int initialLoci = int(loci_ind(0));
  int nCells = HabitatMap.ncol() * HabitatMap.nrow();
  
  int n = int_1;
  
  //// now get within while loop
  int iterations = int_100; //int_1e6;
  IntegerVector events; //define high because used everywhere; 
  
  ////////////////////////////////////////////////////////////////////////////////////
  // while loop
  
  while((loci_ind.size() < terrSizeMax) & ((n - 1) < iterations)){ 
    
    ///////////////////////////////////////////////
    //find potential cells
    
    //loci_x and loci_y, before when was based on XY coord
    //Rcout << "Rcout loci_ind start while: " << std::endl << loci_ind << std::endl;
    List potentials = UniqFreeAdjCellsRandOrd(loci_ind, availCellsMat, 1, 1);
    // sanity check
    //return(potentials);
    IntegerVector potentials_AdjX = potentials["AdjX"];
    IntegerVector potentials_AdjY = potentials["AdjY"];
    IntegerVector potentials_CellNum = potentials["CellNum"];
    //Rcout << "Rcout potentials_CellNum: " << std::endl << potentials_CellNum << std::endl;
    // potential cells kept
    int nPot = potentials_AdjX.size();
    
    if(nPot > 0){
      // decided to put everything in a single loop
      for(int p = 0; p<nPot; p++){
        if(loci_ind.size() < terrSizeMax){
          loci_ind.push_back(potentials_CellNum(p));
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
  
  Rcout << "Rcout loci_ind" << std::endl << loci_ind << std::endl;
  if(cpp_as_output(0) == 0){
    //initialLoci = initialLoci + 1;
    for(int i = 0; loci_ind.size(); i++){
      loci_ind(i) = loci_ind(i) + 1;
    }
  }
  
  List L_return_spreadCpp = List::create(Named("initialLoci") = initialLoci,
                                         _["CellNum"] = loci_ind 
  );
  
  return L_return_spreadCpp;
  
}// end of function


////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
// building a map of all potential territories from sequence of position
List TerrMapping(
    IntegerMatrix HabitatMap, // different from availCellUpdateRas ?
    IntegerMatrix availCellsMat,
    int terrSizeMax,
    int terrSizeMin,
    IntegerVector terrCentreCellNum,
    IntegerMatrix terrMap,
    LogicalVector cpp_as_input,
    LogicalVector cpp_as_output
){
  
  if(cpp_as_input(0) == 0){
    for(int i = 0; i<terrCentreCellNum.size(); i++){
      terrCentreCellNum(i) = terrCentreCellNum(i) - 1;
    }
  }
  
  int int_0 = 0;
  
  List TerrCentreFullCoords = CellNumtoRowCol(terrCentreCellNum, availCellsMat, 1, 1);
  IntegerVector TerrCentreFullX = TerrCentreFullCoords["x_coords"];
  IntegerVector TerrCentreFullY = TerrCentreFullCoords["y_coords"];
  //Rcout << "TerrCentreFullCoords" << std::endl << TerrCentreFullCoords << std::endl;
  
  for(int c = 0; c<terrCentreCellNum.size(); c++){
    
    //Rcout << "cell treated " << std::endl << c << std::endl;
    
    int Xcurrent = TerrCentreFullX(c);
    int Ycurrent = TerrCentreFullY(c);
    
    List potTerrSpread = spreadGB(availCellsMat,
                                  Ycurrent,
                                  Xcurrent,
                                  terrSizeMax,
                                  1, 
                                  1);
    
    IntegerVector potTerrCellNum = potTerrSpread["CellNum"];
    List potTerrFullCoords = CellNumtoRowCol(potTerrCellNum, availCellsMat, 1, 1);
    IntegerVector potTerrX = potTerrFullCoords["x_coords"];
    IntegerVector potTerrY = potTerrFullCoords["y_coords"];
    
    //Rcout << "potTerrCellNum" << std::endl << potTerrCellNum << std::endl;
    
    //Rcout << "Spread done " << std::endl << c << std::endl;
    
    if(potTerrCellNum.size() >= terrSizeMin){ // if at least 9km2
      //Rcout << "Terr big enough " << std::endl << c << std::endl;
      for(int n = 0; n < potTerrCellNum.size(); n++){
        terrMap(c, n) = potTerrCellNum(n);
        availCellsMat(potTerrY(n), potTerrX(n)) = 0;
      }
    }
    
    //Rcout << "Terr rec done " << std::endl << c << std::endl;
    
  }  
  
  List L_return_terrMap = List::create(Named("terrMap") = terrMap,
                                       _["availCellsMat"] = availCellsMat 
  );
  
  return L_return_terrMap;
  
}
