
////////////////////////////////////////////////////////////////////////////////
  
  // [[Rcpp::export]]
// building a map of all potential territories from sequence of position
IntegerMatrix TerrMapping(
  IntegerMatrix HabitatMap, // different from availCellUpdateRas ?
    IntegerMatrix availCellsUpdatedRas,
  int terrSizeMax,
  int terrSizeMin,
  IntegerVector terrCentreCellNum,
  IntegerMatrix terrMap
){
  
  List TerrCentreFullCoords = CellNumtoRowCol(terrCentreCellNum, availCellsUpdatedRas);
  
  for(int c = 0; c<terrCentreCellNum.size(); c++){
    
    List potTerrCellNum = spreadGB(availCellsUpdatedRas, //IntegerMatrix availCellsUpdatedRas,
                                   TerrCentreFullCoords["y_coords"](c), //int DispFem_lastDispY,
                                   TerrCentreFullCoords["x_coords"](c), //int DispFem_lastDispX,
                                   terrSizeMax); //int terrSize) 

if(potTerrCellNum["CellNum"].size() >= terrSizeMin){ // if at least 9km2
  for(int n = 0; n < potTerrCellNum["CellNum"].size(); n++){
    list potTerrCellsFullCoord =  CellNumtoRowCol(cellNum = potTerrCellNum["CellNum"](n), Matrix = availCellsUpdatedRas);
    TerrMap(c, n) = potTerrCellsFullCoord["CellNum"](n);
  }
}

}  
  
  List L_return_terrMap = List::create(Named("terrMap") = terrMap,
                                       _["availCellsUpdatedRas"] = availCellsUpdatedRas 
  );
  
  return L_return_terrMap;
  
  }