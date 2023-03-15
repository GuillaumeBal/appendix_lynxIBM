
reqdPkgs = list("NetLogoR", "testthat", "SpaDES", "raster", "randomcoloR", "data.table", "dplyr", "doBy")
reqdPkgs %>% unlist %>% sapply(., FUN = function(x){require(x, character.only = TRUE)})
require(magrittr)

#############################
## Analyze the simulations ##
#############################
# Read all simulated files

pathFiles <- "C:/Users/gbal/Desktop/lynx.ibm/appendix_lynxIBM/module/calibration/calibr1" %T>% dir
listSim <- list.files(pathFiles)[1:10]#[2:201] # remove the first file (cache)
nSim <- length(listSim)
lastYear <- 50

load(paste0(pathFiles, "/", listSim[3]))

sim <- lynxIBMrun
sim$._knownObjects
sim$._simPrevs
sim %>% names
structure(sim$par)
P(sim, "paths")

