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
// cell number to coords in cpp
List CellNumtoRowCol(IntegerVector cellNums, IntegerMatrix Matrix){
  IntegerVector rowNum(cellNums.size());
  IntegerVector colNum(cellNums.size());
  int nRowMat = Matrix.nrow();
  int nColMat = Matrix.ncol();
  for(int i = 0; i<cellNums.size(); i++){
    rowNum(i) = trunc(cellNums(i) / nColMat);
    colNum(i) = cellNums(i) - rowNum(i) * nColMat - 1;
  }
  List coordsRCNum = List::create(Named("x_coords") = colNum,
                                  _["y_coords"] = rowNum,
                                  _["cellNums"] = cellNums);
  return coordsRCNum;
}

// [[Rcpp::export]]
// adjacent cells coordinates
List UniqFreeAdjCellsRandOrd(IntegerVector cellNum, IntegerMatrix Matrix){
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
  // change cell num into XY
  List all_coords = CellNumtoRowCol(cellNum = cellNum, Matrix = Matrix);
  IntegerVector x_coords = all_coords["x_coords"];
  IntegerVector y_coords = all_coords["y_coords"];
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
  List L_return = List::create(Named("AdjX") = AdjX_left,
                               _["AdjY"] = AdjY_left,
                               _["CellInd"] = CellInd_left);
  // return L_return;
  return L_return;
}
