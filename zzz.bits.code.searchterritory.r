terrSizeName <- 
  of(world = sim$popDist, 
     agents = patchHere(world = sim$popDist, # donnes les dernière coords
                        turtles = turtle(turtles = sim$lynx, who = searchingFemID)))

dispFem[1 , ]
patchHere(world = sim$popDist, 
          turtles = turtle(turtles = sim$lynx, who = searchingFemID))


terrDT <- spread(landscape = availCellsUpdatedRas,
                 loci = cellFromPxcorPycor(world = sim$habitatMap,
                                           pxcor = searchingFemCell[, 1],
                                           pycor = searchingFemCell[, 2]),
                 spreadProb = availCellsUpdatedRas, maxSize = terrSize, returnIndices = TRUE,
                 quick = TRUE)

terrCells <- unique(terrDT$indices)


?cellFromPxcorPycor
?spread


body(spread)
