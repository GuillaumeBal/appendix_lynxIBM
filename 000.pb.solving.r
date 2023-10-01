require(raster)
avail.mat.check <- avail.mat.check <- read.table('gb.avail.matrix.txt') %>% as.matrix
crash.point <- read.table('crash.index.txt') %>% unlist
#raster(avail.mat.check %>% as.matrix()) %>% plot

sim$availCellsRas %>% plot()

cell.crash.cpp <- cell.index[crash.point %>% unlist()]
coords.crash <- CellNumtoRowCol(cellNum = cell.crash.cpp, Matrix = avail.mat.check,
                                cpp_as_input =  1, cpp_as_output = 1)
coords.crash
points(coords.crash$x_coords + 1, coords.crash$y + 1, pch = 16, col = 'red')

UniqFreeAdjCellsRandOrd(
  cellNum = cell.crash.cpp, Matrix = avail.mat.check,
  cpp_as_input = 1, cpp_as_output = 0
)

spreadGB.outputs <-
  spreadGB(
    availCellsMat = sim$habitatMap %>% as.matrix(),#avail.mat.check,#IntegerMatrix availCellsMat,
    YCoordinate = coords.crash$y_coords + 1,
    XCoordinate = coords.crash$x_coords + 1,
    terrSizeMax = terrSize.max,
    cpp_as_input = 0,
    cpp_as_output = 1
  )
spreadGB.outputs$CellNum %>% length

check.spread <- 
  spread(landscape = avail.mat.check,# sim$availCellsRas,
         loci = cellFromPxcorPycor(world = sim$habitatMap,
                                   pxcor = coords.crash$x_coords + 1,
                                   pycor = coords.crash$y_coords + 1),
         spreadProb = sim$availCellsRas, #avail.mat.check,# 
         maxSize = terrSize.max, returnIndices = TRUE,
         quick = TRUE)
check.spread %>% dim()

check.spread.xy <- terra::xyFromCell(cell = check.spread$indices, object = sim$availCellsRas)

plot(sim$availCellsRas,
     xlim = c(coords.crash$x_coords + 1 - 7, coords.crash$x_coords + 1 + 7),
     ylim = c(coords.crash$y_coords + 1 - 7, coords.crash$y_coords + 1 + 7))
points(check.spread.xy,  pch = 16, col = 'blue', cex = 2)
points(coords.crash$x_coord +1 , coords.crash$y +1, pch = 16, col = 'red', cex = 2)
