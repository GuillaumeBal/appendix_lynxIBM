count.dim <- 1
dim.allcells.1 <- 1
while(dim.allcells.1 <= 96){
  x.picked <- sample.int(sim$habitatMap@maxPxcor, 1)
  y.picked <- sample.int(sim$habitatMap@maxPycor, 1)
  allCells <-  spread.short.gb(
    landscape = availCellsUpdatedRas,
    loci = cellFromPxcorPycor(world = sim$habitatMap,
                              pxcor = x.picked,
                              pycor = y.picked),
    spreadProb = availCellsUpdatedRas, 
    maxSize = terrSize,
    returnIndices = TRUE,
    quick = TRUE)
  dim.allcells.1 <- dim(allCells)[1]
  print(paste(count.dim, dim.allcells.1, sep = ' / '))
  count.dim <- count.dim + 1
  if(count.dim>100) break()
}
#allCells
