rm(list = ls())

timeframe = as.POSIXlt(c(NA, NA))
timeunit = "year"
citation = list("citation.bib")
documentation = deparse(list("README.txt", "lynxIBM.Rmd"))
reqdPkgs = list("NetLogoR", "testthat", "SpaDES", "raster",
                "randomcoloR", "data.table", "dplyr", "doBy")

require(magrittr)

"C:/Users/gbal/Desktop/lynx.ibm/appendix_lynxIBM/" %>% setwd()

reqdPkgs %>% unlist %>% sapply(., FUN = function(x){require(x, character.only = TRUE)})

?of
body(of)
?standardGeneric
