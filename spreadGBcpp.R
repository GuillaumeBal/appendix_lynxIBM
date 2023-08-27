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

sourceCpp("spreadGBcpp.cpp")

x.picked.base <- 50#sample.int(sim$habitatMap@maxPxcor, 1)
y.picked.base <- 380#sample.int(sim$habitatMap@maxPycor, 1)
x.picked <- x.picked.base - 1 # cpp correction as starts at 0
y.picked <- y.picked.base - 1 

rm('outputs.cpp')
outputs.cpp <- try(
  spreadGB(
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
    deadLynxColl = sim$deadLynxColl[[time(sim, "year")[1]]][,] %>% `[`( ,c('who', "steps", "heading", "lastDispX", "lastDispY")),
    deadDisp = sim$deadDisp[,],
    TerrMap = sim$terrMap@.Data[,],  # bits for searchterritory,
    availCellsUpdatedRas = sim$availCellsRas %>% as.matrix,
    popDist = sim$popDist@.Data[,],
    coreTerrSizeFAlps = coreTerrSizeFAlps,
    coreTerrSizeFJura = coreTerrSizeFJura,
    coreTerrSizeFVosgesPalatinate = coreTerrSizeFVosgesPalatinate,
    coreTerrSizeFBlackForest = coreTerrSizeFBlackForest,
    returnDistances = FALSE,
    allowOverlap = FALSE,
    DispFem_lastDispY = y.picked,# DispFem_lastDispY(f) in full cpp 
    DispFem_lastDispX = x.picked,# DispFem_lastDispX(f) in full cpp
    terrSize = terrSize # in full cpp
  )
)

outputs.cpp

rm('outputs.loop')
n_try <- 10000
record.freq <- 5
outputs.loop <- list()

rec <- 1
for(i in 1:n_try){ # run several times to check for potential indexing issues that are sometimes fine a few times
  x.picked <- sample.int(sim$habitatMap@maxPxcor, 1) - 1
  y.picked <- sample.int(sim$habitatMap@maxPycor, 1) - 1 
  if(i == 1){ 
    y.x.piched.df <- data.frame("y.picked" = y.picked, "x.picked" = x.picked)
    write.table(
      y.x.piched.df,
      file = 'x.y.picked.txt')
  }
  #y.x.piched.df <- rbind(y.x.piched.df, c(y.picked, x.picked))
  write.table(c(y.picked, x.picked), file = 'x.y.picked.txt', append = TRUE)
  outputs.loop.current <- 
    #tryCatch(
    try(
      spreadGB(
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
        allowOverlap = FALSE,
        DispFem_lastDispY = y.picked,# DispFem_lastDispY(f) in full cpp 
        DispFem_lastDispX = x.picked,# DispFem_lastDispX(f) in full cpp
        terrSize = terrSize # in full cpp
      )#,
      #error = function(e) 
    )
  if(i %in% seq(0, n_try, record.freq)){
    outputs.loop[[rec]] <- outputs.loop.current
    outputs.loop[[rec]]$x.picked = x.picked
    outputs.loop[[rec]]$y.picked = y.picked
    rec <- rec + 1
  }
  if(i %in% seq(0, n_try, n_try/100)) print(i)
  #print(outputs.loop[[i]])
  #if(outputs.loop$MatInd %>% length %>% `==`(0)) stop()
}

