require(inline)
require(Rcpp)
require(RcppArmadillo)
require(tryCatchLog)
#require(futile.logger)
#options("tryCatchLog.write.error.dump.file" = TRUE)

# work on search territory part ====================================================

lynx.gb <- sim$lynx[,]
lynx.gb$steps <- NULL # drop it because created within cpp
habitatMap.gb <- sim$habitatMap@.Data[,] # x coords are on columns and y on lines
terrMap.gb <- sim$terrMap@.Data[,]
availCellsUpdatedRas.gb <- sim$availCellsRas %>% as.matrix() # in fact before update
popDist.gb <- sim$popDist@.Data
trick <- c(1, 1)

sourceCpp("dispersalGBcpp.cpp")

outputs.cpp <- try(
  dispersalGB(
    lynx_r = lynx.gb,
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
    deadLynxColl = sim$deadLynxColl[[time(sim, "year")[1]]][,] %>% `[`(,c('who', "steps", "heading", "lastDispX", "lastDispY")),
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
# outputs.cpp
# outputs.cpp %>% sapply(., length)
# outputs.cpp$lynx_dispNow %>% length
# outputs.cpp$Lynx_status %>% length
# outputs.cpp$Lynx_lastDispX %>% length


outputs.loop <- list()
for(i in 1:500){ # run several times to check for potential indexing issues that are sometimes fine a few times
  outputs.loop[[i]] <- 
    #tryCatch(
    try(
      dispersalGB(
        lynx_r = lynx.gb,#[lynx.gb$who != 1868, ],
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
        deadLynxColl = sim$deadLynxColl[[time(sim, "year")[1]]][,] %>% `[`(,c('who', "steps", "heading", "lastDispX", "lastDispY")), # did not bother to do the others
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
      )#,
      #error = function(e) e
    )
  #print(i)
  #print(outputs.loop[[i]])
  #if(outputs.loop$MatInd %>% length %>% `==`(0)) stop()
}

#outputs.loop
#outputs.loop[[1]] %>% sapply(., FUN = length)

# check for run return before issue data =======================================
error.cases <- 
  outputs.loop %>% 
  sapply(.,
         function(x){
           #names(x)[1] == 'int_100'
           x[1] %>% as.character(.) %>% substr(., 1, 5)
         }) %>% `==`('Error') %>% which 

error.cases %>% length 

okay.cases <- 
  outputs.loop %>% 
  sapply(.,
         function(x){
           #names(x)[1] == 'int_100'
           x[1] %>% as.character(.) %>% substr(., 1, 5)
         }) %>% `!=`('Error') %>% which 

matching.size.issue <- outputs.loop %>% 
  sapply(.,
         function(x){
           #if(x[1] %>% `==` ('mat_chosen')){ #substr(., 1, 5) %>%
             if(x[1] %>% substr(., 1, 5) %>% `!=` ('Error')){ #substr(., 1, 5) %>%
             #if(length(x$ChosenCell_who) != x$nDispLeft){
             if(length(unique(c(x$ChosenCellsYesCorr_who, x$ChosenCellsNoCorr_who))) != x$nDispLeft){
             #if(length(x$dispersers_who) != length(x$Mat_Chosen)){
               return(TRUE)
             }
           }
           return(FALSE)
         }) %>% which 
matching.size.issue

outputs.loop[[matching.size.issue[1]]]

