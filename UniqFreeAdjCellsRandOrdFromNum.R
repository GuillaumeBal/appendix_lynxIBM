require(inline)
require(Rcpp)
require(RcppArmadillo)
require(tryCatchLog)
#require(futile.logger)
#options("tryCatchLog.write.error.dump.file" = TRUE)

# work on search territory part ====================================================

lynx.gb <- sim$lynx[,]
lynx.gb$steps <- NULL # drop it because created within cpp
habitatMap.gb <- sim$habitatMap@.Data[,] # x coords are on columns and y on lines
terrMap.gb <- sim$terrMap@.Data[,]
availCellsUpdatedRas.gb <- sim$availCellsRas %>% as.matrix() # in fact before update
popDist.gb <- sim$popDist@.Data
trick <- c(1, 1)

sourceCpp("UniqFreeAdjCellsRandOrdFromNum.cpp")

x.picked.base <- 47
y.picked.base <- 101

x.picked <- x.picked.base - 1 #<- 428
y.picked <- y.picked.base - 1  #<- 523

nColMat <-ncol(sim$availCellsRas %>% as.matrix) 
nRowMat <- nrow(sim$availCellsRas %>% as.matrix)

cellnum.picked <- nColMat * (nRowMat - (y.picked.base)) + (x.picked.base);

rm('outputs.cpp')
outputs.cpp <- try(
  UniqFreeAdjCellsRandOrd(
    Matrix = sim$availCellsRas %>% as.matrix,
    cellNum = c(181515, 181942, 182372, 181944, 182373, 182371)#cellnum.picked
    #y_coords = y.picked,# DispFem_lastDispY(f) in full cpp 
    #x_coords = x.picked
  )
)

outputs.cpp

rm('outputs.loop')
n_try <- 20000
record.freq <- 10
outputs.loop <- list()

rec <- 1
file.remove('x.y.picked.txt')
for(i in 1:n_try){ # run several times to check for potential indexing issues that are sometimes fine a few times
  x.picked <- sample.int(sim$habitatMap@maxPxcor, 1) - 1
  y.picked <- sample.int(sim$habitatMap@maxPycor, 1) - 1
  if(i == 1){ 
    y.x.piched.df <- data.frame("y.picked" = y.picked, "x.picked" = x.picked)
    write.table(
      y.x.piched.df,
      file = 'x.y.picked.txt')
  }
  #y.x.piched.df <- rbind(y.x.piched.df, c(y.picked, x.picked))
  write.table(c(y.picked, x.picked), file = 'x.y.picked.txt', append = TRUE)
  outputs.loop.current <- 
    #tryCatch(
    try(
      UniqFreeAdjCellsRandOrd(
        Matrix = sim$availCellsRas %>% as.matrix,
        y_coords = y.picked,# DispFem_lastDispY(f) in full cpp 
        x_coords = x.picked
      )#,
      #error = function(e) 
    )
  if(i %in% seq(0, n_try, record.freq)){
    outputs.loop[[rec]] <- outputs.loop.current
    outputs.loop[[rec]]$x.picked = x.picked
    outputs.loop[[rec]]$y.picked = y.picked
    rec <- rec + 1
    saveRDS(outputs.loop, file = '0.outputs.to.check.RDS')
  }
  if(i %in% seq(0, n_try, n_try/100)) print(i)
  #print(outputs.loop[[i]])
  #if(outputs.loop$MatInd %>% length %>% `==`(0)) stop()
}

