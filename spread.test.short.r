require(Rcpp)
require(inline)
sourceCpp("spreadGBcpp.cpp")

x.picked.cpp <- x.picked -1
y.picked.cpp <- y.picked -1

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
  DispFem_lastDispY = y.picked.cpp,# DispFem_lastDispY(f) in full cpp 
  DispFem_lastDispX = x.picked.cpp,# DispFem_lastDispX(f) in full cpp
  terrSize = terrSize # in full cpp
)
