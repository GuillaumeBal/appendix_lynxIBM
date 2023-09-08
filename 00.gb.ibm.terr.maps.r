require(inline)
require(Rcpp)
require(RcppArmadillo)
require(tryCatchLog)

# need a cpp to deal with all combination

map.avail.matrix <- sim$availCellsRas %>% as.matrix
cell.map.max <- dim(map.avail.matrix)[1] * dim(map.avail.matrix)[2]
cell.start <- runif(1, 0, cell.map.max) %>% round 
cell.index <- c(seq(cell.start, cell.map.max, by = 3),
                seq(0, cell.start, by = 3))
terrSize.max <- 97
terrMatrix <- matrix(-100, nrow = cell.index %>% length, ncol = terrSize.max)

sourceCpp("00.terrMapping.cpp")

outputs.spread <- try(
  spreadGB(
    availCellsMat = sim$availCellsRas %>% as.matrix,
    XCoordinate = 131 - 1,#cell.map.max - 1,
    YCoordinate = 380 - 1,#terrSize.max, # in full cpp
    terrSizeMax = terrSize.max
  )
)

time.start <- Sys.time()
outputs.terr <- try(
  TerrMapping(
    HabitatMap = sim$availCellsRas %>% as.matrix, #different from availCellUpdateRas ?
    availCellsMat = sim$availCellsRas %>% as.matrix,
    terrSizeMax = 97,
    terrSizeMin = 9 ,
    terrCentreCellNum = cell.index[1:500] - 1,
    terrMap = terrMatrix
  )
)
difftime(time.start, Sys.time())

outputs.terr$terrMap %>% table()
