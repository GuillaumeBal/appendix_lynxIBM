require(inline)
require(Rcpp)
require(RcppArmadillo)

# work on search territory part ====================================================

disp.gb <- sim$lynx[sim$lynx$status == 'disp'][ , ]
habitatMap.gb <- sim$habitatMap@.Data[,]
terrMap.gb <- sim$terrMap@.Data[,]
availCellsUpdatedRas.gb <- sim$availCellsRas %>% as.matrix() # in fact before update
popDist.gb <- sim$popDist@.Data
trick <- c(1,1)
source("dispersalGBcppfunction.r")
dispersalGB(
  disperser = disperser[,],
  sMaxPs = sMaxPs,
  trick = trick)

terrMap.gb %>% table(., useNA = 'always')
sim$terrMap@.Data[,] %>% table
availCellsUpdatedRas.gb %>% class
popDist.gb %>% class
sim$habitatMap@.Data %>% c %>% table
dispFem
popDist.gb[dispFem[1, "lastDispX"], dispFem[1, "lastDispY"]]
popDist.gb[297, 425]
popDist.gb %>% dim


