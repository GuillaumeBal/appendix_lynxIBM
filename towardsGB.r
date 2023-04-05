require(inline)
require(Rcpp)
require(RcppArmadillo)

disp.gb.agentmat <- sim$lynx[sim$lynx$status == 'disp']

disp.gb.agentmat %>% colnames()

# work on data frame format ====================================================

disp.gb <- disp.gb.agentmat[ , ]
habitatMap.gb <- sim$habitatMap@.Data[,]
terrMap.gb <- sim$terrMap@.Data[ , ]
sourceCpp("towardsGB.cpp")
towardsGB(Disp = disp.gb, 
          HabitatMap = habitatMap.gb, 
          TerrMap = terrMap.gb)
