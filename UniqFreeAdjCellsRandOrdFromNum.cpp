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
  //Rcout << "Rcout cells with duplicates" << std::endl << CellInd << std::endl;
  // find unique index and then pick just one of each
  IntegerVector UniqCell = unique(CellInd);
  //Rcout << "Rcout UniqCell" << std::endl << UniqCell << std::endl;
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
