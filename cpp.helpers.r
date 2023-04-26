require(Rcpp)
sourceCpp("cpp.helpers.cpp") 

# count number times cells vector match character string
N_Eq_Str(letters, 'a')

# give back order of vector of integer
x <- sample.int(10, 13, replace = TRUE)
order(x)
IntOrderIndex(x)  + 1 # need +1 because c++ starts at 0

# return index cells with values above criteria
y <- sample.int(5, 13, replace = TRUE)
y
index.to.keep <- WhichAbove(y, 3) + 1 #because c++ stats at 0
index.to.keep 
IntVecSubIndex(y, WhichAbove(y, 3)) # here no +1 as back in cpp

# return index num cells whom values are withing a set
z <- sample.int(10, 13, replace = TRUE)
z.sub <- z[sample.int(length(z), 4, replace = FALSE)] %>% unique
which(z %in% z.sub)
WhichInSetInt(z, z.sub) + 1

# return a random line for index from Integer vactor with rep
zz <- sample.int(10, 20, replace = TRUE)
zz
unique(zz) %>% length
zz[IntPosOneOfEach(zz) + 1] #%>% length

# check Whixh Equal for 1
OneAndZero <- rbinom(13, 1, .5)
OneAndZero
WhichEqual(OneAndZero, 1) +1 

# check towards
x_to = c(296,297,298,296,298,296,297,298,297)
y_to = c(426,426,426,425,425,424,424,424,425)
x_cur = 297
y_cur = 425
towards_simple(x_to = x_to, y_to = y_to, x_cur = x_cur, y_cur = y_cur)

# change heading value
changeHeading(90)

# ==============================================================================

sourceCpp("cpp.helpers.cpp")
#return unique adjacent cells 
my.mat <- matrix(rep(0, 81), nrow = 9) %>% as.matrix()
y_coords = c(4, 4) 
x_coords = c(1, 2)
my.mat[y_coords, x_coords] <- 1
adj_cells <- UniqAdjCells(y_coords = y_coords- 1, x_coords = x_coords- 1, Matrix = my.mat)
adj_cells$AdjX <- adj_cells$AdjX + 1
adj_cells$AdjY <- adj_cells$AdjY + 1 
my.mat[adj_cells$AdjY, adj_cells$AdjX] <- my.mat[adj_cells$AdjY, adj_cells$AdjX] + 1

plot(raster(my.mat))

# ==============================================================================
vec <- 1:10
ShortenIntVec(vec, 8)
length(vec)
