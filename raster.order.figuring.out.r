require(raster)
require(magrittr)
require(NetLogoR)

my.mat <- matrix(0, nrow = 10, ncol = 10)
my.rast <- raster(my.mat)

# regular raster count
my.rast[1] <- 1
my.rast[2] <- 2
my.coords <- c(2,3)
my.rast[my.coords <- c(2,3)] <- 3

# lynx map count ===============================================================
my.rast[cellFromPxcorPycor(world = createWorld(minPxcor = 1, maxPxcor = 10,
                                               minPycor = 1, maxPycor = 10),
                           pxcor = 3,#10,#searchingFemCell[, 1],
                           pycor = 2)] <- 4

# compute same cell number =====================================================
nrow(my.mat) * ncol(my.mat) - 
  ((ncol(my.mat)) * (my.coords[1]) - my.coords[2])

ncol(my.mat) * (nrow(my.mat) - my.coords[1]) + my.coords[2]


plot(my.rast)

# regular row and cell
row.y <- 83 %>% `/`(ncol(my.mat)) %>% round()
col.x <- 83 %>% `-`(row.y * ncol(my.mat) + 1)
row.y
col.x


# from coordinates matrix regular
coords.y <- nrow(my.mat) - 1 - (83  %>% `/`(ncol(my.mat)) %>% round())
coords.x <- 83 %>% `-`(row.y * ncol(my.mat) + 1)
coords.y
coords.x

# from map cell number to map coordinates
y.map <- 1
x.map <- 2
cell.map.coord <- (y.map) * ncol(my.mat) + x.map + 1

# from ncell map
coords.y <- 13 %>% `/`(ncol(my.mat)) %>% ceiling() %>% `-`(1)
coords.x <- 13 %>% `-`(coords.y * ncol(my.mat)) -1
coords.y
coords.x
