require(inline)
require(Rcpp)
require(RcppArmadillo)

disp.gb.agentmat <- sim$lynx[sim$lynx$status == 'disp']

# work on data frame format ====================================================

disp.gb <- disp.gb.agentmat[ , ]
sourceCpp("towardsGB.cpp")
towardsGB(disp.gb)
