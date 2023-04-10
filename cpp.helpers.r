require(Rcpp)
sourceCpp("cpp.helpers.cpp") 

# count number times cells vector match character string
N_Eq_Str(letters, 'a')

# give back order of vector of integer
x <- sample.int(10, 13, replace = TRUE)
order(x)
IntOrderIndex(x)  + 1 # need +1 because c++ starts at 0

