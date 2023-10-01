# cell to check as cpp format
cell.to.check <- sample.int(n = length(cell.index), size = 1000, replace = FALSE)

# gb xy from cpp
gb.cells.coord <- CellNumtoRowCol(cellNum = cell.to.check, Matrix = avail.mat.check,
                                  cpp_as_input = TRUE, cpp_as_output = 0) %>% as.data.frame()
# terra xy from cpp + 1
terra.cells.coord <- terra::xyFromCell(cell = cell.to.check + 1, object = sim$availCellsRas) %>% 
  as.data.frame()

checking.xy.df <- data.frame(
  gb = paste(gb.cells.coord$x_coords, gb.cells.coord$y_coords, sep = '_'), 
  terra = paste(terra.cells.coord$x, terra.cells.coord$y, sep = '_')
)

checking.xy.df %>% head(5)

checking.xy.df$gb %>% 
  `==` (checking.xy.df$terra) %>% table

######################################################################################
# compare the SPREAD functions

cell.spread <- 8268#sample.int(n = length(cell.index), size = 1, replace = FALSE)
cell.spread.xy <- CellNumtoRowCol(cellNum = cell.spread, Matrix = avail.mat.check,
                                  cpp_as_input = TRUE, cpp_as_output = TRUE)

rast.plot.xlim <- c(cell.spread.xy$x_coords + 1 - 10, cell.spread.xy$x_coords + 1 + 10)
rast.plot.ylim <- c(cell.spread.xy$y_coords + 1 - 10, cell.spread.xy$y_coords + 1 + 10)

first_uniq <- UniqFreeAdjCellsRandOrd(
  cellNum = cell.spread, Matrix = sim$availCellsRas %>% as.matrix,#sim$habitatMap %>% as.matrix,
  cpp_as_input = 1, cpp_as_output = 0
)
second_uniq <- UniqFreeAdjCellsRandOrd(
  cellNum = first_uniq$CellNum, Matrix = sim$availCellsRas %>% as.matrix,#sim$habitatMap %>% as.matrix,
  cpp_as_input = 0, cpp_as_output = 0
)
second_uniq$CellNum %>% length
second_uniq$CellNum %>% unique %>% length
third_uniq <- UniqFreeAdjCellsRandOrd(
  cellNum = second_uniq$CellNum, Matrix = sim$availCellsRas %>% as.matrix,#sim$habitatMap %>% as.matrix,
  cpp_as_input = 0, cpp_as_output = 0
)
third_uniq$CellNum %>% length
third_uniq$CellNum %>% unique %>% length
fourth_uniq <- UniqFreeAdjCellsRandOrd(
  cellNum = third_uniq$CellNum, Matrix = sim$availCellsRas %>% as.matrix,#sim$habitatMap %>% as.matrix,
  cpp_as_input = 0, cpp_as_output = 0
)
fourth_uniq$CellNum %>% length
fourth_uniq$CellNum %>% unique %>% length

plot(sim$availCellsRas, xlim = rast.plot.xlim, ylim = rast.plot.ylim)
points(first_uniq$AdjX, first_uniq$AdjY,  pch = 16, col = 'brown', cex = 2)
points(cell.spread.xy$x_coord, cell.spread.xy$y_coords, pch = 16, col = 'red', cex = 2)

spreadGB.outputs <-
  spreadGB(
    availCellsMat =  sim$availCellsRas %>% as.matrix,#sim$habitatMap@.Data,#avail.mat.check,
    YCoordinate = cell.spread.xy$y_coords,
    XCoordinate = cell.spread.xy$x_coords,
    terrSizeMax = terrSize.max,
    cpp_as_input = 1,
    cpp_as_output = 1
  )
spreadGB.outputs %<>% cbind(., CellNumtoRowCol(cellNum = spreadGB.outputs$CellNum, Matrix = avail.mat.check,
                                               cpp_as_input = TRUE, cpp_as_output = 1) %>% as.data.frame())

check.spread <- 
  spread(landscape = avail.mat.check,# sim$availCellsRas,
         loci = cellFromPxcorPycor(world = sim$habitatMap,
                                   pxcor = cell.spread.xy$x_coords + 1,
                                   pycor = cell.spread.xy$y_coords + 1),
         spreadProb = sim$availCellsRas, #avail.mat.check,# 
         maxSize = terrSize.max, returnIndices = TRUE,
         quick = TRUE)

spreadGB.outputs$CellNum %>% unique %>% length
check.spread %>% dim()

check.spread.xy <- terra::xyFromCell(cell = check.spread$indices, object = sim$availCellsRas)

plot(sim$availCellsRas, xlim = rast.plot.xlim, ylim = rast.plot.ylim)
points(check.spread.xy,  pch = 16, col = 'blue', cex = 2)
points(spreadGB.outputs$x_coords, spreadGB.outputs$y_coords,  pch = 16, col = 'brown', cex = 2)
points(cell.spread.xy$x_coord + 1 , cell.spread.xy$y + 1, pch = 16, col = 'red', cex = 2)


