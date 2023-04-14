#microbenchmark(base_R = dist_haversine(xlon, xlat, ylon, ylat),
#               Rcpp = dist_haversine_rcpp(xlon, xlat, ylon, ylat),
#               times = 1000L)