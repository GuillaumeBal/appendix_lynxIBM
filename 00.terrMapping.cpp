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
// adjacent cells coordinates
List UniqFreeAdjCellsRandOrd(IntegerVector cellNum, IntegerMatrix Matrix){
  IntegerVector deltaX = {-1, 0, 1};
  IntegerVector deltaY = {-1, 0, 1};
  int nColMat = Matrix.ncol();
  int nRowMat = Matrix.nrow();
  int int_1 = 1;
  int int_2 = 2;
  IntegerVector AdjX(0);
  IntegerVector AdjY(0);
  IntegerVector CellNum(0);
  int new_y;
  int new_x;
  int new_index;
  // change cell num into XY
  List L_1 = List::create(Named("int_1") = int_1);
  //return(L_1);
  //Rcout << "Rcout before cellnum " << std::endl << int_1 << std::endl;
  List all_coords = CellNumtoRowCol(cellNum = cellNum, Matrix = Matrix);
  List L_2 = List::create(Named("all_coords") = all_coords);
  //return(L_2);
  IntegerVector x_coords = all_coords["x_coords"];
  IntegerVector y_coords = all_coords["y_coords"];
  //Rcout << "Rcout 1 inits" << std::endl << CellNum << std::endl;
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
        bool new_index_exist = IntInSet(new_index, cellNum); // check does not already exist
        if(
          (new_index>=0) &
            (new_y>=0) & (new_y<nRowMat) & (new_x>=0) & (new_x<nColMat) & 
            (new_index_exist == false)){ // add to keep only unoccupied cells
          if((Matrix(nRowMat - 1 - new_y, new_x) == int_1)){
            AdjX.push_back(new_x);
            AdjY.push_back(new_y);
            CellNum.push_back(new_index);
          }
        }
      }
    }
  }
  //Rcout << "Rcout cells with duplicates" << std::endl << CellNum << std::endl;
  // find unique index and then pick just one of each
  IntegerVector UniqCell = unique(CellNum);
  //Rcout << "Rcout UniqCell" << std::endl << UniqCell << std::endl;
  int nUniqCell = UniqCell.size();
  IntegerVector keptCells(nUniqCell);
  for(int i = 0; i<nUniqCell; i++){
    for(int j = 0; j<CellNum.size(); j++){
      if(CellNum(j) == UniqCell(i)){
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
  IntegerVector CellNum_left = CellNum[keptCellsRandOrder];
  //Rcout << "Rcout CellNum_left" << std::endl << CellNum_left << std::endl;
  List L_return = List::create(Named("CellNum") = CellNum_left,
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
    int terrSizeMax
){
  
  ///////////////////////////////////////////////////////////////////////////////////////
  // TRICKS AND FUNCTIONS FOR CODE
  
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
  loci_ind(0) = HabitatMap.ncol() * (HabitatMap.nrow() - (YCoordinate + 1)) + (XCoordinate + 1);
  IntegerVector loci_y(1); loci_y(0) = YCoordinate;
  IntegerVector loci_x(1); loci_x(0) =XCoordinate;
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
    List potentials = UniqFreeAdjCellsRandOrd(loci_ind, availCellsMat);
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
  
  //Rcout << "Rcout n" << std::endl << n << std::endl;
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
    IntegerMatrix terrMap
){
  
  int int_0 = 0;
  
  List TerrCentreFullCoords = CellNumtoRowCol(terrCentreCellNum, availCellsMat);
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
                                  terrSizeMax);
    
    IntegerVector potTerrCellNum = potTerrSpread["CellNum"];
    List potTerrFullCoords = CellNumtoRowCol(potTerrCellNum, availCellsMat);
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
