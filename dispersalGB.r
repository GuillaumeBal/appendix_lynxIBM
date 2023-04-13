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
  pCorr = pCorr,
  nMatMax = nMatMax,
  connectivityMap = sim$connectivityMap@.Data[,],
  roadMortMap = sim$roadMortMap@.Data[,],
  corrFactorDisp = corrFactorDisp,
  floorTimeSim = floor(time(sim))[1],
  startSimYear = start(sim, "year")[1],
  ncoll_ncoll = sim$nColl$ncoll,
  ncoll_time = sim$nColl$time,
  deadLynxColl = sim$deadLynxColl[[time(sim, "year")[1]]][,]
)
outputs.cpp

for(i in 1:500){ # run several times to check for potential indexing issues that are sometimes fine a few times
  outputs.loop <- try(
    dispersalGB(
      lynx = lynx.gb,
      sMaxPs = sMaxPs,
      HabitatMap = habitatMap.gb,
      pMat = pMat, #round(1/9, 2)
      pCorr = pCorr,
      nMatMax = nMatMax,
      connectivityMap = sim$connectivityMap@.Data[,],
      roadMortMap = sim$roadMortMap@.Data[,],
      corrFactorDisp = corrFactorDisp,
      floorTimeSim = floor(time(sim))[1],
      startSimYear = start(sim, "year")[1],
      ncoll_ncoll = sim$nColl$ncoll,
      ncoll_time = sim$nColl$time,
      deadLynxColl = sim$deadLynxColl[[time(sim, "year")[1]]][,]
    ) 
  )
  print(outputs.loop)
  #if(outputs.loop$MatInd %>% length %>% `==`(0)) stop()
}



