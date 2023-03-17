rm(list = ls())

library(igraph)
library(SpaDES.core)
library(raster)

# SORTIES SAUVEES DONC PLUS BESOIN DE FAIRE TOURNER

# Where is located your appendix_lynxIBM repo
pathRepo <- "C:/Users/gbal/Desktop/lynx.ibm/" # change to your path

# stuff changed by GB
lynx.data <- 
  "appendix_lynxIBM/module/inputs/listLynxInitPopSub.RData"
  #"appendix_lynxIBM/module/inputs/listLynxInitPop.RData"
n.years.run <- 1 #50

# Define the paths
moduleDir <- file.path(paste0(pathRepo, "appendix_lynxIBM/module"))
inputDir <- file.path(paste0(pathRepo, moduleDir, "inputs")) %>% 
  reproducible::checkPath(create = TRUE)
outputDir <- file.path(paste0(pathRepo, moduleDir, "outputs"))
cacheDir <- file.path(paste0(pathRepo, outputDir, "cache"))

# Module parameters
times <- list(start = 1, end = n.years.run + 1) #need to add 1 to have desired number years
parameters <- list(
  .plotInitialTime = NA, # no plotting
  .plotInterval = NA,
  testON = TRUE # activate the inside-function tests 
  )

# Input files
# Habitat quality raster
habMapSpaDES <- raster(paste0(pathRepo, "appendix_lynxIBM/module/inputs/habMap.tif"))
# Collision probabilities raster
collProbSpaDES <- raster(paste0(pathRepo, "appendix_lynxIBM/module/inputs/collProb.tif"))
# List of 500 different initial population (SpatialPointsDataFrame)
load(paste0(pathRepo, lynx.data))
popInitSpaDES <- listLynxInitPop[[sample(1:length(listLynxInitPop), 1)]] # sample one initial population
# Population areas raster
fourPopSpaDES <- raster(paste0(pathRepo, "appendix_lynxIBM/module/inputs/fourPop.tif"))

modules <- list("lynxIBM")
objects = c("habMapSpaDES", "collProbSpaDES", "popInitSpaDES", "fourPopSpaDES")
paths <- list(
  cachePath = cacheDir,
  modulePath = moduleDir,
  inputPath = inputDir,
  outputPath = outputDir
)

# Initialize the module
lynxIBMinit <- simInit(times = times, params = list(lynxIBM = parameters), modules = modules,
                objects = objects, paths = paths)

# Run one simulation
run.start <- Sys.time()
lynxIBMrun <- spades(lynxIBMinit, debug = TRUE)
run.end <- Sys.time()

difftime(run.start, run.end)

P(lynxIBMinit)
lynx.sim.gb <- list('lynxIBMrun' = lynxIBMrun, 'lynxIBMinit' = lynxIBMinit)
save(lynx.sim.gb, file = paste0(pathRepo, 'appendix_lynx/IBMlynx.sim.gb.Rdata'))
# Save the simulation output
#save(mySimOut, file = paste0(pathRepo, "appendix_lynxIBM/module/outputs/lynxIBMrun_", sample(1:100000,1), ".RData"))
