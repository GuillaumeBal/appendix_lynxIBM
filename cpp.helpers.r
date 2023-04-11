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



