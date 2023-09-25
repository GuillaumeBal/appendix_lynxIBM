require(raster)
avail.mat.check <- avail.mat.check <- read.table('gb.avail.matrix.txt') %>% as.matrix
crash.point <- read.table('crash.index.txt') %>% unlist
#raster(avail.mat.check %>% as.matrix()) %>% plot

sim$availCellsRas %>% plot()

cell.crash.cpp <- cell.index[crash.point %>% unlist() %>% `+`(0)]
coords.crash <- CellNumtoRowCol(cellNum = cell.crash.cpp, Matrix = avail.mat.check)
coords.crash
points(coords.crash$x_coords, coords.crash$y, pch = 16, col = 'red')
avail.mat.check %>% dim

UniqFreeAdjCellsRandOrd(
  cellNum = cell.crash.cpp, Matrix = avail.mat.check
)

spreadGB.outputs <-
spreadGB(
  availCellsMat = avail.mat.check,#IntegerMatrix availCellsMat,
  YCoordinate = coords.crash$y_coords,
  XCoordinate = coords.crash$x_coords,
  terrSizeMax = terrSize.max
)


CellNumtoRowCol(spreadGB.outputs$CellNum, Matrix = avail.mat.check) # some pb of negative x_coords that should not happen


points(spreadGB.outputs$)


spread(landscape = sim$availCellsUpdatedRas,
       loci = cellFromPxcorPycor(world = sim$habitatMap,
                                 pxcor = coords.crash$x_coords + 1,
                                 pycor = coords.crash$y_coords + 1),
       spreadProb = sim$availCellsUpdatedRas, maxSize = terrSize.max, returnIndices = TRUE,
       quick = TRUE)



#terr.updated <- 
sapply(outputs.terr, FUN = 
         function(x){
           x$terrMap %>% `==`(-100) %>% table
         }) %>% do.call(., what = cbind)
outputs.terr[[3]]$terrMap[, 1] %>% `==`(-100) %>% table


225225
197337
