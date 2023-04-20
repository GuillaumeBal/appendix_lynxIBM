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
  lynx_list = lynx.gb %>% as.data.frame %>% as.list,
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
  deadLynxColl = sim$deadLynxColl[[time(sim, "year")[1]]][,],
  deadDisp = sim$deadDisp[,],
  TerrMap = sim$terrMap@.Data[,],  # bits for searchterritory,
  availCellsUpdatedRas = sim$availCellsRas %>% as.matrix,
  popDist = sim$popDist@.Data[,],
  coreTerrSizeFAlps = coreTerrSizeFAlps,
  coreTerrSizeFJura = coreTerrSizeFJura,
  coreTerrSizeFVosgesPalatinate = coreTerrSizeFVosgesPalatinate,
  coreTerrSizeFBlackForest = coreTerrSizeFBlackForest,
  returnDistances = FALSE,
  allowOverlap = FALSE
)
outputs.cpp

outputs.loop <- list()
for(i in 1:500){ # run several times to check for potential indexing issues that are sometimes fine a few times
  outputs.loop[[i]] <- 
    try(
      dispersalGB(
        lynx = lynx.gb,
        lynx_list = lynx.gb %>% as.data.frame %>% as.list,
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
        deadLynxColl = sim$deadLynxColl[[time(sim, "year")[1]]][,],
        deadDisp = sim$deadDisp[,],
        TerrMap = sim$terrMap@.Data[,],  # bits for searchterritory,
        availCellsUpdatedRas = sim$availCellsRas %>% as.matrix,
        popDist = sim$popDist@.Data[,],
        coreTerrSizeFAlps = coreTerrSizeFAlps,
        coreTerrSizeFJura = coreTerrSizeFJura,
        coreTerrSizeFVosgesPalatinate = coreTerrSizeFVosgesPalatinate,
        coreTerrSizeFBlackForest = coreTerrSizeFBlackForest,
        returnDistances = FALSE,
        allowOverlap = FALSE
      ) 
    )
  #print(i)
  #print(outputs.loop[[i]])
  #if(outputs.loop$MatInd %>% length %>% `==`(0)) stop()
}

outputs.loop[[1]]

outputs.loop %>% 
  sapply(., 
         function(x){
           sum(x$deathRoad)
         }) %>%  table

outputs.loop %>% 
  sapply(., 
         function(x){
           lynx_new.length <- x['lynx_xcor_new'] %>% length
           lynx.length <-  x['lynx_xcor']  %>% length
           return(lynx_new.length < lynx.length)
         }) %>% table
lynx.gb$maleID