require(inline)
require(Rcpp)
require(RcppArmadillo)

# work on search territory part ====================================================

lynx.gb <- sim$lynx[,]
lynx.gb$steps <- NULL # drop it because created within cpp
habitatMap.gb <- sim$habitatMap@.Data[,] # x coords are on columns and y on lines
terrMap.gb <- sim$terrMap@.Data[,]
availCellsUpdatedRas.gb <- sim$availCellsRas %>% as.matrix() # in fact before update
popDist.gb <- sim$popDist@.Data
trick <- c(1, 1)
sourceCpp("dispersalGBcpp.cpp")
outputs.cpp <- dispersalGB(
  lynx = lynx.gb,
  sMaxPs = sMaxPs,
  HabitatMap = habitatMap.gb,
  pMat = pMat, #round(1/9, 2)  #pMat
  pCorr = pCorr
)
outputs.cpp
for(i in 1:500){ # run several times to check for potential indexing issues that are sometimes fine a few times
  try(print(
    dispersalGB(
      lynx = lynx.gb,
      sMaxPs = sMaxPs,
      HabitatMap = habitatMap.gb,
      pMat = pMat, #round(1/9, 2)
      pCorr = pCorr
    ) 
  ))
}
