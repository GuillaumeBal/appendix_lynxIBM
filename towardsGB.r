require(inline)
require(Rcpp)
require(RcppArmadillo)

disp.gb <- sim$lynx[sim$lynx$status == 'disp']
disp.gb <- disp.gb[,]

sourceCpp("towardsGB.cpp")

towardsGB(disp.gb)


disp.gb[,] %>% class
