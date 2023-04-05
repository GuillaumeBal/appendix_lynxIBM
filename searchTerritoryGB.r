require(inline)
require(Rcpp)
require(RcppArmadillo)

# work on search territory part ====================================================

disp.gb <- sim$lynx[sim$lynx$status == 'disp'][ , ]
habitatMap.gb <- sim$habitatMap@.Data[,]
terrMap.gb <- sim$terrMap@.Data[,]
availCellsUpdatedRas.gb <- sim$availCellsRas %>% as.matrix() # in fact before update
popDist.gb <- sim$popDist@.Data
sourceCpp("searchTerritoryGB.cpp")
searchTerritoryGB(Disp = disp.gb, 
          HabitatMap = habitatMap.gb, 
          TerrMap = terrMap.gb,
          availCellsUpdatedRas = availCellsUpdatedRas.gb,
          popDist = popDist.gb,
          coreTerrSizeFAlps = coreTerrSizeFAlps,
          coreTerrSizeFJura = coreTerrSizeFJura,
          coreTerrSizeFVosgesPalatinate = coreTerrSizeFVosgesPalatinate,
          coreTerrSizeFBlackForest = coreTerrSizeFBlackForest)


terrMap.gb %>% table(., useNA = 'always')
sim$terrMap@.Data[,] %>% table
availCellsUpdatedRas.gb %>% class
popDist.gb %>% class
sim$habitatMap@.Data %>% c %>% table
dispFem
popDist.gb[dispFem[1, "lastDispX"], dispFem[1, "lastDispY"]]
popDist.gb[297, 425]
popDist.gb %>% dim


