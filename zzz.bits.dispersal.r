x <- sample.int(10, 25, replace = TRUE)
xsort <- sort(unique(x))
Index <- rep(NA, length(x))
p = 1
for(i in 1:length(xsort)){
  for(j in 1:length(x)){
    if(x[j] == xsort[i]){
      Index[p] = j
      p = p+1;
    }
  }
}

x[Index]
sum(x[Index]) == sum(x)