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
outputs.cpp

outputs.loop <- list()
for(i in 1:500){ # run several times to check for potential indexing issues that are sometimes fine a few times
  outputs.loop[[i]] <- 
    #tryCatch(
    try(
      dispersalGB(
        lynx = lynx.gb,#[lynx.gb$who != 1868, ],
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

outputs.loop

outputs.loop[[1]] %>% sapply(., FUN = length)

# check for run return before issue data =======================================
no.corres.cases <- 
  outputs.loop %>% 
  sapply(.,
         function(x){
           #names(x)[1] == 'int_100'
           x[1] %>% as.character(.) %>% substr(., 1, 5)
         }) %>% `==`('Error') %>% which 

outputs.loop %>% sapply(., function(x){
  sum(x$lynx_who_new == 0) == x$nDisp_new
}) %>% table

outputs.loop %>% sapply(., function(x){
  sum(x$deathRoad)
}) %>% table


no.corres.cases
picked.to.look <- 4995
output.looked <- outputs.loop[[picked.to.look]]
output.looked
output.looked$lynx_who %>% `==` (output.looked$who_issue) %>% which
output.looked$dispersers_who %in% output.looked$lynx_who

lynx.gb[lynx.gb$who == 1868, ]

# number of steps 
outputs.loop %>% 
  sapply(.,
         function(x){
           x$step
         }) %>% `[`(no.corr.cases) %>%  table


#==============================================================
outputs.loop %>% 
  sapply(.,
         function(x){
           x$deadLynxColl %>% class
         }) %>% `==`("list") %>% which


outputs.loop %>% 
  sapply(.,
         function(x){
           length(x$lynx_who) %>% table %>% table
         }) %>%  table

outputs.loop %>% 
  sapply(.,
         function(x){
           length(x$dispersers_who_new) %>% table %>% table
         }) %>%  table


# outputs.loop %>% 
#     sapply(.,
#            function(x){
#              substr(x$message, nchar(x$message)-20+1, nchar(x$message))
#            }) %>%  table