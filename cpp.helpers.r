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
