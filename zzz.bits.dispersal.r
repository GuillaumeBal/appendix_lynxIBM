outputs.cpp


nextCellsType_x <- rep(NA, outputs.cpp$nCellsDispLeft)
nextCellsType_y <- rep(NA, outputs.cpp$nCellsDispLeft)
nextCellsType_ind <- rep(NA, outputs.cpp$nCellsDispLeft)
nextCellsType_hab <- rep(NA, outputs.cpp$nCellsDispLeft)
p = 0
for(i in 1:length(outputs.cpp$CellsDisp_ind)){
  if((outputs.cpp$Mat_chosen[outputs.cpp$CellsDisp_ind[i]+1] == 1) & (outputs.cpp$CellsDisp_hab[i] == 2)){
    nextCellsType_x[p] = outputs.cpp$CellsDisp_x[i]
    nextCellsType_y[p] = outputs.cpp$CellsDisp_y[i]
    nextCellsType_ind[p] = outputs.cpp$CellsDisp_ind[i]
    nextCellsType_hab[p] = outputs.cpp$CellsDisp_hab[i]
    p = p + 1
  }
  if((outputs.cpp$Mat_chosen[outputs.cpp$CellsDisp_ind[i]+1] == 0) &
     (( outputs.cpp$CellsDisp_hab[i] == 3) | ( outputs.cpp$CellsDisp_hab[i] == 4))){
    nextCellsType_x[p] = outputs.cpp$CellsDisp_x[i]
    nextCellsType_y[p] = outputs.cpp$CellsDisp_y[i]
    nextCellsType_ind[p] = outputs.cpp$CellsDisp_ind[i]
    nextCellsType_hab[p] = outputs.cpp$CellsDisp_hab[i]
    p = p + 1
  }
}
