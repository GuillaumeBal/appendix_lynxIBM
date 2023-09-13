require(inline)
require(Rcpp)
require(RcppArmadillo)
require(tryCatchLog)

# need a cpp to deal with all combination

map.avail.matrix <- sim$availCellsRas %>% as.matrix
update.cells <- 'no' # no to be able to pick saved ones to find issues within the script
if(update.cells == 'yes'){
  terrSize.max <- 97
  cell.map.max <- dim(map.avail.matrix)[1] * dim(map.avail.matrix)[2]
  cell.start <- runif(1, 0, cell.map.max) %>% round 
  cell.index <- c(seq(cell.start, cell.map.max, by = 3),
                  seq(0, cell.start, by = 3))
  terrMatrix <- matrix(-100, nrow = cell.index %>% length, ncol = terrSize.max)
  #saved details so that i can reload them
  cells.details = list(
    'cell.map.max' = cell.map.max,
    'cell.start' = cell.start,
    'cell.index' = cell.index,
    'terrMatrix' = terrMatrix,
    'terrSize.max' = terrSize.max
  )
  saveRDS(cells.details, 'cells.details.RDS')
}else{
  attach(readRDS('cells.details.RDS'))
}
terrSize.min <-  9

sourceCpp("00.terrMapping.cpp")

#time.start <- Sys.time()
# system.time(
#   outputs.terr <- try(
#     TerrMapping(
#       HabitatMap = sim$availCellsRas %>% as.matrix, #different from availCellUpdateRas ?
#       availCellsMat = sim$availCellsRas %>% as.matrix,
#       terrSizeMax = terrSize.max,
#       terrSizeMin = terrSize.min,
#       terrCentreCellNum = cell.index - 1,#cell.index - 1,#[1:100] - 1,
#       terrMap = terrMatrix
#     )
#   )
# )
#difftime(time.start, Sys.time())

gb.avail.matrix <- sim$availCellsRas %>% as.matrix

#recording before crash
outputs.terr <- list()
for(i in 1:length(cell.index)){
  write.table(i, file = 'crash.index.txt')
  outputs.terr[[i]] <- try(
    TerrMapping(
      HabitatMap = gb.avail.matrix, #different from availCellUpdateRas ?
      availCellsMat = gb.avail.matrix,
      terrSizeMax = terrSize.max,
      terrSizeMin = terrSize.min,
      terrCentreCellNum = cell.index[i] - 1,#cell.index - 1,#[1:100] - 1,
      terrMap = terrMatrix
    )
  )
  #print(outputs.terr)
  if(length(outputs.terr[[i]]) == 2){
  gb.avail.matrix <- outputs.terr[[i]]$availCellsMat
  write.table(gb.avail.matrix, file = 'gb.avail.matrix.txt')
  }else{
    break()
  }
}

outputs.terr %>% sapply(., length) %>% table
outputs.terr[[1]]$availCellsMat

#outputs.terr$terrMap %>% `[`(1: 100, 1:10)
#outputs.terr$terrMap %>% table() %>% table()

# old check of spread function outputs
# outputs.spread <- try(
#   spreadGB(
#     availCellsMat = sim$availCellsRas %>% as.matrix,
#     XCoordinate = 131 - 1,#cell.map.max - 1,
#     YCoordinate = 380 - 1,#terrSize.max, # in full cpp
#     terrSizeMax = terrSize.max
#   )
# )