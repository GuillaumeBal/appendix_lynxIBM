rm(list = ls())

require(magrittr)

wd <- "C:/Users/gbal/Desktop/appendix_lynxIBM/module/inputs" %T>% setwd() 
dir()

load("listLynxInitPop.RData")

listLynxInitPop[[10]]
listLynxInitPop[[10]] %T>% class() %>% dim() 

listLynxInitPopSub <- 
  sapply(listLynxInitPop, 
         function(x){
           kept <- sample.int(n = x %>% dim(.) %>% `[`(1),
                              size = x %>% dim(.) %>% `[`(1) %>% `/`(., 10) %>% round())
           x <- x[kept , ]
           return(x)
         })

# need to keep same names list in IBM
listLynxInitPop <- listLynxInitPopSub 
ls()
rm(list = c("wd", "listLynxInitPopSub"))

save(listLynxInitPop, file = "listLynxInitPopSub.RData")


