my.mat <- matrix(0, nrow = 10, ncol = 10)
my.rast <- raster(my.mat)

# regular raster count
my.rast[1] <- 1
my.rast[2] <- 2
my.coords <- c(2,3)
my.rast[my.coords <- c(2,3)] <- 3

# lynx map count
my.rast[cellFromPxcorPycor(world = createWorld(minPxcor = 1, maxPxcor = 10,
                                              minPycor = 1, maxPycor = 10),
                           pxcor = 3,#10,#searchingFemCell[, 1],
                           pycor = 2)] <- 4

nrow(my.mat) * ncol(my.mat) - 
  ((ncol(my.mat)) * (my.coords[1]) - my.coords[2])

ncol(my.mat) * (nrow(my.mat) - my.coords[1]) + my.coords[2]


plot(my.rast)




