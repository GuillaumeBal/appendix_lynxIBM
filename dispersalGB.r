require(inline)
require(Rcpp)
require(RcppArmadillo)

# work on search territory part ====================================================

lynx.gb <- sim$lynx[,]
lynx.gb$steps <- NULL # drop it because created within cpp
habitatMap.gb <- sim$habitatMap@.Data[,]
terrMap.gb <- sim$terrMap@.Data[,]
availCellsUpdatedRas.gb <- sim$availCellsRas %>% as.matrix() # in fact before update
popDist.gb <- sim$popDist@.Data
trick <- c(1, 1)
sourceCpp("dispersalGBcpp.cpp")
dispersalGB(
  lynx = lynx.gb,
  sMaxPs = sMaxPs,
  HabitatMap = habitatMap.gb,
  pMat = round(1/9, 2)  #pMat
)
for(i in 1:500){
  print(
    dispersalGB(
      lynx = lynx.gb,
      sMaxPs = sMaxPs,
      HabitatMap = habitatMap.gb,
      pMat = round(1/9, 2)  #pMat
    ) 
  )
}
