# figuring out coordinates for loci ################################################## 

loci = as.integer(cellFromPxcorPycor(world = sim$habitatMap,
                                     pxcor = searchingFemCell[, 1],
                                     pycor = searchingFemCell[, 2]))


which(sim$habitatMap@pCoords[,"pxcor"] == searchingFemCell[, 1] &
        sim$habitatMap@pCoords[,"pycor"] == searchingFemCell[, 2])

loci
searchingFemCell[, 1]
searchingFemCell[, 2]
sim$habitatMap %>% dim
dim.hab <- dim(sim$habitatMap)
sim$habitatMap@pCoords

matrix(c(1,2,3,4), ncol = 2) %>% c

# pxcor
dim.hab[2] * (dim.hab[1] - searchingFemCell[, 2]) +  searchingFemCell[, 1]


####################################################################################
dim.1.spredf <- 1
while(dim.1.spredf <= 1){
x.picked <- sample.int(sim$habitatMap@maxPxcor, 1)
y.picked <- sample.int(sim$habitatMap@maxPycor, 1)
spread.res <- spread(landscape = availCellsUpdatedRas,
       loci = cellFromPxcorPycor(world = sim$habitatMap,
                                 pxcor = x.picked,#10,#searchingFemCell[, 1],
                                 pycor = y.picked),#520),#searchingFemCell[, 2]),
       spreadProb = availCellsUpdatedRas, maxSize = terrSize, returnIndices = TRUE,
       quick = TRUE)
dim.1.spredf <- spread.res %>% dim %>% `[`(1)
}

######################################################################################
# check spread outputs 

# first things, car def ========================================
outputs.cpp$loci_ind
loci

outputs.cpp$spreadsDT_spreads %>% length
spreadsDT %>% nrow

outputs.cpp$spreadsDT_spreads %>% `==`(spreadsDT$spreads) %>% table
outputs.cpp$spreadsDT_spreads[loci]
outputs.cpp$spreadsDT_spreads %>% sum
spreadsDT$spreads[loci]

#
