require(raster)
avail.mat.check <- avail.mat.check <- read.table('gb.avail.matrix.txt') %>% as.matrix
crash.point <- read.table('crash.index.txt') %>% unlist
#raster(avail.mat.check %>% as.matrix()) %>% plot

sim$availCellsRas %>% plot()

cell.crash.cpp <- cell.index[crash.point %>% unlist() %>% `+`(0)]
coords.crash <- CellNumtoRowCol(cellNum = cell.crash.cpp, Matrix = avail.mat.check,
                                cpp_as_input =  1, cpp_as_output = 0)
coords.crash
points(coords.crash$x_coords, coords.crash$y, pch = 16, col = 'red')

terra::xyFromCell(cell = cell.crash.cpp + 1, object = sim$availCellsRas)


UniqFreeAdjCellsRandOrd(
  cellNum = cell.crash.cpp, Matrix = avail.mat.check,
  cpp_as_input = 1, cpp_as_output = 0
)

spreadGB.outputs <-
  spreadGB(
    availCellsMat = avail.mat.check,#IntegerMatrix availCellsMat,
    YCoordinate = coords.crash$y_coords,
    XCoordinate = coords.crash$x_coords,
    terrSizeMax = terrSize.max,
    cpp_as_input = 0,
    cpp_as_output = 1
  )
spreadGB.outputs$CellNum %>% length











check.spread.coord <- CellNumtoRowCol(spreadGB.outputs$CellNum, Matrix = avail.mat.check) # some pb of negative x_coords that should not happen
check.spread.coord %<>% do.call(., what = 'rbind')
check.spread.coord[, check.spread.coord['x_coords', ] < 0]
check.spread.coord
dim(avail.mat.check)

nRowMat <- nrow(avail.mat.check)
nColMat <- ncol(avail.mat.check)
cellNum <- check.spread.coord[, check.spread.coord['x_coords', ] < 0][1]
rowNum = (nRowMat - (max(0, floor((cellNum  + 1)) / nColMat)) + 1) - 1
colNum = ((cellNum + 1) - (nRowMat - (rowNum + 1)) * nColMat) - 1 

points(spreadGB.outputs$)

check.spread <- 
  spread(landscape = avail.mat.check,# sim$availCellsRas,
         loci = cellFromPxcorPycor(world = sim$habitatMap,
                                   pxcor = coords.crash$x_coords + 1,
                                   pycor = coords.crash$y_coords + 1),
         spreadProb = avail.mat.check,#sim$availCellsUpdatedRas, 
         maxSize = terrSize.max, returnIndices = TRUE,
         quick = TRUE)
check.spread %>% dim()






nColMat <- ncol(avail.mat.check)
nRowMat <- nrow(avail.mat.check)
base_rowNum =  ceiling((cell.crash.cpp + 1) / nColMat);
rowNum_xy = (nRowMat - base_rowNum + 1) - 1 ; 
colNum_xy = ((cell.crash.cpp + 1) - ((base_rowNum - 1) * nColMat)) - 1
rowNum_xy
colNum_xy

#terr.updated <- 
sapply(outputs.terr, FUN = 
         function(x){
           x$terrMap %>% `==`(-100) %>% table
         }) %>% do.call(., what = cbind)
outputs.terr[[3]]$terrMap[, 1] %>% `==`(-100) %>% table


225225
197337
