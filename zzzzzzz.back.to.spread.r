cells <- c(182370,182805,183662,182376,182799,184092,184091,183664,181520,181519,182374,183229,184090,182377,181946,179799,182802,182804,181943,181082,181090,180656,181512,181942,181518,181517,182801,183227,181948,183228,183235,182800,183233,181944,182372,182797,181089,183236,181949,182373,184093,182375,181085,181945,183665,181941,184089,182371,181515,183661,181083,183656,181516,183663,180227,182806)

cells %>% length

cells %>% duplicated() %>% table()

outputs.cpp$events %>% unique %>% length

outputs.cpp$events


x.picked.base <- sample.int(sim$habitatMap@maxPxcor, 1)
y.picked.base <- sample.int(sim$habitatMap@maxPycor, 1)

while(dim(spread.res)[1]< 2 | dim(spread.res)[1] ==  terrSize){
  
  #x.picked.base <- sample.int(sim$habitatMap@maxPxcor, 1)
  #y.picked.base <- sample.int(sim$habitatMap@maxPycor, 1)
  
  spread.res <- spread(landscape = availCellsUpdatedRas,
                       loci = cellFromPxcorPycor(world = sim$habitatMap,
                                                 pxcor = x.picked.base,#10,#searchingFemCell[, 1],
                                                 pycor = y.picked.base),#520),#searchingFemCell[, 2]),
                       spreadProb = availCellsUpdatedRas,
                       maxSize = terrSize, 
                       returnIndices = TRUE,
                       quick = TRUE)
}

# plot female potential terr zone
availCellsUpdatedRas.sub <- crop(x = availCellsUpdatedRas, y = extent(x.picked.base - 13, x.picked.base+13, y.picked.base - 13, y.picked.base+13))
plot(availCellsUpdatedRas.sub)

# points sarah's points for terr
spread.cells <- CellNumtoRowCol(cellNum = spread.res$indices, Matrix = availCellsUpdatedRas.gb)
points(spread.cells$x_coords + 1, spread.cells$y_coords + 1,  col = 'orange', lwd = 2)

# points mine
spread.cells.gb <- CellNumtoRowCol(cellNum = outputs.cpp$loci_ind, Matrix = availCellsUpdatedRas.gb)
points(spread.cells.gb$x_coords + 1, spread.cells.gb$y_coords + 1,  col = 'purple', lwd = 2)

# put female current loc
points(x.picked.base, y.picked.base, col = 'red', lwd = 2) 

