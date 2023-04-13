rm(list = ls())

require(magrittr)
wd <- "C:/Users/gbal/Desktop/lynx.ibm/appendix_lynxIBM" %T>% setwd()
# agent matrix
load("lynx.sim.gb.Rdata")
#lynx.sim.gb$lynxIBMrun
#attach(lynx.sim.gb)

load('../calibr.sarah/calibr1/lynxIBMrun_10319.RData')
lynxIBMinit <- lynx.sim.gb$lynxIBMinit

# IBM linked packages
reqdPkgs = list("NetLogoR", "testthat", "SpaDES", "raster", "randomcoloR", "data.table", "dplyr", "doBy")
reqdPkgs %>% unlist %>% sapply(., FUN = function(x){require(x, character.only = TRUE)})

# some data from outputs
sim <- lynxIBMrun#lynxIBMinit

# arguments Sarah to get back from inits part of loaded run data
sMaxPs <- P(lynxIBMinit)$lynxIBM$sMaxPs #45
testON <- P(lynxIBMinit)$lynxIBM$testON #TRUE
pMat <- P(lynxIBMinit)$lynxIBM$pMat #0.03
pCorr <- P(lynxIBMinit)$lynxIBM$pCorr #0.5
nMatMax <- P(lynxIBMinit)$lynxIBM$nMatMax #10
corrFactorDisp <- P(lynxIBMinit)$lynxIBM$corrFactorDisp #250 
# within searchTerritory
coreTerrSizeFAlps <- P(lynxIBMinit)$lynxIBM$coreTerrSizeFAlps #97
coreTerrSizeFJura <- P(lynxIBMinit)$lynxIBM$coreTerrSizeFJura #126
coreTerrSizeFVosgesPalatinate <- P(lynxIBMinit)$lynxIBM$coreTerrSizeFVosgesPalatinate #126
coreTerrSizeFBlackForest <- P(lynxIBMinit)$lynxIBM$coreTerrSizeFBlackForest #126
terrSizeMAlps <- P(lynxIBMinit)$lynxIBM$terrSizeMAlps #159
terrSizeMJura <- P(lynxIBMinit)$lynxIBM$terrSizeMJura #270
terrSizeMVosgesPalatinate <- P(lynxIBMinit)$lynxIBM$terrSizeMVosgesPalatinate #270
terrSizeMBlackForest <- P(lynxIBMinit)$lynxIBM$terrSizeMBlackForest #270

#"C:/Users/gbal/Desktop/lynx.ibm/searchterritory.func.r" %>% 
#  source()

# here is where the function start , dispersal <- function(sim) { ===================================================
# stops before executing part to code within cpp

disperser <- NLwith(agents = sim$lynx, var = "status", val = "disp")
nDisp <- NLcount(disperser)
nonDisperser <- other(agents = sim$lynx, except = disperser)
nonDisperserID <- nonDisperser@.Data[, "who"]
sim$deadDisp <- rbind(sim$deadDisp, data.frame(nDisp = nDisp, nDispDeadColl = 0, nDispDeadDaily = 0, time = floor(time(sim))[1]))

#if(nDisp != 0) {
  
  # Number of steps dispersers can walk during the day
  steps <- sample(x = 1:sMaxPs, size = nDisp, prob = sim$Ps, replace = TRUE) # GB_change P(sim)$sMaxPs
  stepsDisp <- 1:max(steps)
  #disperser <- NLset(turtles = disperser, agents = disperser, var = "steps", val = steps)
  disperser@.Data[, "steps"] <- as.numeric(steps) # faster
  
  #stop()
  
  #{for(step in stepsDisp) {
  step <- 2  
  #if(NLcount(disperser) != 0){ # update of the disperser at each "step" loop
  
  dispersingIndNMatMax <- noTurtles()
  
  # Dispersers which have not reach their max number of steps yet
  infoDisp <- of(agents = disperser, var = c("who", "status", "steps"))
  dispersingID <- infoDisp[step <= infoDisp[, "steps"] & infoDisp[, "status"] == "disp", "who"]
  dispersingInd <- turtle(turtles = disperser, who = dispersingID)
  nonDispersingInd <- other(agents = disperser, except = dispersingInd)
  
  # Test
  if(testON == TRUE) { #GB_change P(sim)$testON
    infoDispersingInd <- of(agents = dispersingInd, var = c("status", "steps"))
    expect_true(all(infoDispersingInd[, "status"] == "disp"))
    expect_true(all(infoDispersingInd[, "steps"] >= step))
    infoNonDispersingInd <- of(agents = nonDispersingInd, var = c("status", "steps"))
    if(NLcount(nonDispersingInd) != 0){
      expect_true(all(infoNonDispersingInd[, "status"] == "res" |
                        infoNonDispersingInd[, "steps"] < step))
    }
  }
  
  #if(NLcount(dispersingInd) != 0){
  
  # Next step in one of the 9 (or less) cells = 8 (or less) neighboring cells + current location
  neighCells <- NetLogoR::neighbors(world = sim$habitatMap, agents = dispersingInd, nNeighbors = 8,
                                    torus = FALSE)
  currCells <- patchHere(world = sim$habitatMap, turtles = dispersingInd)
  currCells <- cbind(currCells, id = 1:NLcount(dispersingInd))
  nextCells <- rbind(neighCells, currCells)
  
  # First, selection of cell type with preference for breeding or dispersal over matrix
  nextCells <- cbind(nextCells,
                     cellType = of(world = sim$habitatMap,
                                   agents = nextCells[, c("pxcor", "pycor")]))
  nextCellsDT <- as.data.table(nextCells)
  cellTypeFreqFull <- nextCellsDT[ , count := .N, by = list(id, cellType)]
  cellTypeFreqFull <- cellTypeFreqFull[, c("id", "cellType", "count")]
  #setkey(cellTypeFreqFull) # not necessary with dplyr::distinct()
  #cellTypeFreq <- as.matrix(unique(cellTypeFreqFull)) # cbind is faster than as.matrix
  #cellTypeFreqUnik <- unique(cellTypeFreqFull)
  cellTypeFreqUnik <- dplyr::distinct(cellTypeFreqFull) # faster than unique
  cellTypeFreq <- cbind(id = cellTypeFreqUnik$id, cellType = cellTypeFreqUnik$cellType,
                        count = cellTypeFreqUnik$count)
  
  # Test
  if(testON == TRUE) { # GB_ change P(sim)$testON
    expect_equivalent(sum(cellTypeFreq[, "count"]), NROW(nextCells))
  }
  
  withMatrix <- unique(cellTypeFreq[cellTypeFreq[, "cellType"] == 2, "id"])
  withoutMatrix <- unique(cellTypeFreq[!cellTypeFreq[, "id"] %in% withMatrix, "id"])
  pChooseMatrix <- pMat * cellTypeFreq[cellTypeFreq[, "id"] %in% withMatrix & # GB_change P(sim)$pMat
                                         cellTypeFreq[, "cellType"] == 2, "count"]
  choseMatrix <- rbinom(n = length(pChooseMatrix), size = 1, prob = pChooseMatrix)
  onlyMat <- unique(cellTypeFreq[, "id"])[
    !unique(cellTypeFreq[, "id"]) %in%
      unique(cellTypeFreq[cellTypeFreq[, "cellType"] %in% c(3, 4) , "id"])]
  nextCellsType <- rbind(nextCells[nextCells[, "cellType"] == 2 &
                                     nextCells[, "id"] %in%
                                     c(withMatrix[choseMatrix == 1], onlyMat), ],
                         nextCells[nextCells[, "cellType"] %in% c(3, 4) &
                                     nextCells[, "id"] %in% c(withoutMatrix,
                                                              withMatrix[choseMatrix == 0]), ])
  
  # Second, selection of cell based on the direction
  if(step == 1) { # no correlation
    
    nextCellsTypeDT <- as.data.table(nextCellsType)
    nextCellsTypeDTsampled <- nextCellsTypeDT[nextCellsTypeDT[, .I[sample(.N,1)], by = id]$V1]
    # chosenCells <- as.matrix(nextCellsTypeDTsampled) # the following line is fater than this one
    chosenCells <- cbind(pxcor = nextCellsTypeDTsampled$pxcor, pycor = nextCellsTypeDTsampled$pycor,
                         id = nextCellsTypeDTsampled$id, cellType = nextCellsTypeDTsampled$cellType)
    chosenCells <- chosenCells[order(chosenCells[, "id"]), , drop = FALSE]
    
    # Test
    if(testON == TRUE) { #P(sim)$testON
      expect_equivalent(NROW(chosenCells), length(unique(chosenCells[, "id"])))
      expect_equivalent(NROW(chosenCells), length(unique(nextCellsType[, "id"])))
    }
   
  # END step == 1   
  } else { # movement correlation
    
    dispCell <- cbind(patchHere(world = sim$habitatMap, turtles = dispersingInd),
                      id = 1:NLcount(dispersingInd))
    colnames(dispCell)[c(1, 2)] <- c("pxcorHere", "pycorHere")
    #stop()
    nextCellsType <- merge(nextCellsType, dispCell)
    nextCellsType <- nextCellsType[order(nextCellsType[, "id"]), , drop = FALSE]
    
    # Individuals with a correlated movement
    probCorr <- rbinom(n = NLcount(dispersingInd), size = 1, prob = pCorr) # GB_change P(sim)$pCorr
    
    # Individuals without a correlated movement
    noCorr <- nextCellsType[nextCellsType[, "id"] %in%
                              unique(nextCellsType[, "id"])[probCorr == 0], ]
    chosenCellsNoCorr <- cbind(pxcor = numeric(), pycor = numeric(),  id = numeric(),
                               cellType = numeric())
    if(length(noCorr) != 0) {
      noCorrDT <- as.data.table(noCorr)
      # chosenCellsNoCorr <- as.matrix(noCorrDT[noCorrDT[, .I[sample(.N,1)], by = id]$V1]) # the following two lines are faster than this single one line
      chosenCellsNoCorrDT <- noCorrDT[noCorrDT[, .I[sample(.N,1)], by = id]$V1]
      chosenCellsNoCorr <- cbind(id = chosenCellsNoCorrDT$id, pxcor = chosenCellsNoCorrDT$pxcor,
                                 pycor = chosenCellsNoCorrDT$pycor, cellType = chosenCellsNoCorrDT$cellType,
                                 pxcorHere = chosenCellsNoCorrDT$pxcorHere, pycorHere = chosenCellsNoCorrDT$pycorHere)
      chosenCellsNoCorr <- chosenCellsNoCorr[, c("pxcor", "pycor", "id", "cellType")]
    }
    
    #stop()
    
    # Individuals with a correlated movement
    yesCorr <- unique(nextCellsType[, "id"])[probCorr == 1]
    chosenCellsYesCorrSelect <- cbind(pxcor = numeric(), pycor = numeric(),  id = numeric(),
                                      cellType = numeric())
    if(length(yesCorr) != 0) {
      nextCellsTypeDir <- cbind(nextCellsType[nextCellsType[, "id"] %in% yesCorr, ],
                                dir = as.numeric(NA))
      for(eachYesCorr in yesCorr) {
        indYessCorr <- turtle(turtles = dispersingInd, who = dispersingID[eachYesCorr])
        #stop()
        #source("b2.towards.defunc.r")
        #start.towards <- Sys.time()
        dirCells <- towards(agents = indYessCorr, agents2 = cbind(
          pxcor = nextCellsTypeDir[nextCellsTypeDir[, "id"] == eachYesCorr, "pxcor"],
          pycor = nextCellsTypeDir[nextCellsTypeDir[, "id"] == eachYesCorr, "pycor"]),
          torus = FALSE)
        #print(paste('towards', start.towards - Sys.time()))
        nextCellsTypeDir[nextCellsTypeDir[, "id"] == eachYesCorr, "dir"] <-
          round(subHeadings(angle1 = indYessCorr@.Data[, "heading"], angle2 = dirCells,
                            range360 = TRUE))
        
      }
      
      # Test
      if(testON == TRUE) { # GB_change P(sim)$testON
        expect_true(all(nextCellsTypeDir[, "dir"] %in% c(0, 45, 90, 135, 180, 225, 270, 315)))
        expect_equivalent(length(yesCorr), length(unique(nextCellsTypeDir$id)))
      }
      
      # Rank the direction as preferences
      nextCellsTypeDir <- cbind(nextCellsTypeDir, prefDir = 1)
      nextCellsTypeDir[nextCellsTypeDir[, "dir"] %in% c(45, 315), "prefDir"] <- 2
      nextCellsTypeDir[nextCellsTypeDir[, "dir"] %in% c(90, 270), "prefDir"] <- 3
      nextCellsTypeDir[nextCellsTypeDir[, "dir"] %in% c(135, 225), "prefDir"] <- 4
      nextCellsTypeDir[nextCellsTypeDir[, "dir"] == 180, "prefDir"] <- 5
      nextCellsTypeDir[nextCellsTypeDir[, "pxcor"] == nextCellsTypeDir[, "pxcorHere"]
                       & nextCellsTypeDir[, "pycor"] == nextCellsTypeDir[, "pycorHere"],
                       "prefDir"] <- 3
      nextCellsTypeDirDT <- as.data.table(nextCellsTypeDir)
      chosenCellsYesCorrDT <- nextCellsTypeDirDT[nextCellsTypeDirDT[, .I[sample(.N,1)],
                                                                    by = c("id","prefDir")]$V1]
      # Select the cell with the smallest rotation
      #chosenCellsYesCorrSelect <- as.matrix(chosenCellsYesCorrDT[, .SD[which.min(prefDir)], by = id]) # the two lines after are faster than these two lines
      #chosenCellsYesCorrSelect <- chosenCellsYesCorrSelect[, c("pxcor", "pycor", "id", "cellType")]
      chosenCellsYesCorrSelectDT <- chosenCellsYesCorrDT[, .SD[which.min(prefDir)], by = id]
      chosenCellsYesCorrSelect <- cbind(pxcor = chosenCellsYesCorrSelectDT$pxcor, pycor = chosenCellsYesCorrSelectDT$pycor,
                                        id = chosenCellsYesCorrSelectDT$id, cellType = chosenCellsYesCorrSelectDT$cellType)
    }
    
    # Regroup the individuals that had not a correlation movement and the ones that had
    chosenCells <- rbind(chosenCellsNoCorr, chosenCellsYesCorrSelect)
    chosenCells <- chosenCells[order(chosenCells[, "id"]), , drop = FALSE]
    
    # Test
    if(testON == TRUE) { # GB_change P(sim)$testON
      expect_equivalent(NROW(chosenCells), length(unique(chosenCells[, "id"])))
      expect_equivalent(NROW(chosenCells), length(unique(nextCellsType[, "id"])))
      expect_equivalent(NROW(chosenCells), NLcount(dispersingInd))
    }
  }
  
  # Lynx memory
  chosenMat <- chosenCells[chosenCells[, "cellType"] == 2, , drop = FALSE]
  if(NROW(chosenMat) != 0) {
    dispersingIndMat <- turtle(turtles = dispersingInd, who = dispersingID[chosenMat[, "id"]])
    dispersingIndMatnSteps <- dispersingIndMat@.Data[, "nMat"]
    #dispersingInd <- NLset(turtles = dispersingInd, agents = dispersingIndMat, var = "nMat",
    #                       val = dispersingIndMatnSteps + 1)
    dispersingInd@.Data[match(dispersingIndMat@.Data[, "who"], 
                              dispersingInd@.Data[, "who"]), "nMat"] <- as.numeric(dispersingIndMatnSteps + 1) # faster
    
    # Use memory to find a dispersal cell
    if(sum(dispersingIndMatnSteps + 1 == nMatMax) != 0) { # GB_change P(sim)$nMatMax
      chosenCells <- as.data.frame(chosenCells)
      chosenCells[chosenCells[, "cellType"] == 2, c("pxcor", "pycor")][
        dispersingIndMatnSteps + 1 == nMatMax, ] <- dispersingIndMat[ # GB_change P(sim)$nMatMax
          dispersingIndMatnSteps + 1 == nMatMax]@.Data[, c("lastDispX", "lastDispY")] # GB_change P(sim)$nMatMax
      
      # Test
      if(testON == TRUE) { # GB_change P(sim)$testON
        expect_equivalent(of(agents = NLwith(agents = dispersingInd, var = "nMat", val = nMatMax), # GB_change P(sim)$nMatMax
                             var = "who"),
                          of(agents = dispersingIndMat[dispersingIndMatnSteps + 1 == nMatMax], # GB_change P(sim)$nMatMax
                             var = "who"))
        expect_equivalent(of(agents = NLwith(agents = dispersingInd, var = "nMat", val = nMatMax), # GB_change P(sim)$nMatMax
                             var = c("lastDispX", "lastDispY")),
                          dispersingIndMat[dispersingIndMatnSteps + 1 == nMatMax]@.Data[ # GB_change P(sim)$nMatMax
                            , c("lastDispX", "lastDispY")])
      }
      
      # Reset nMat
      dispersingIndNMatMax <- NLwith(agents = dispersingInd, var = "nMat", val = nMatMax) # GB_change P(sim)$nMatMax
      #dispersingInd <- NLset(turtles = dispersingInd, agents = dispersingIndNMatMax, var = "nMat",
      #                       val = 0)
      dispersingInd@.Data[match(dispersingIndNMatMax@.Data[, "who"], dispersingInd@.Data[, "who"]), "nMat"] <- as.numeric(0) # faster
      
      # Test
      if(testON == TRUE) { # GB_change P(sim)$testON
        lastDispCell <- of(agents = NLwith(agents = dispersingInd, var = "nMat", val = nMatMax), # GB_change P(sim)$nMatMax
                           var = c("lastDispX", "lastDispY"))
        expect_equivalent(sum(is.na(lastDispCell)), 0)
        expect_true(all(of(agents = dispersingInd, var = "nMat") < nMatMax)) # GB_change P(sim)$nMatMax
      }
    }
  }
  
  # Update "lastDispX" and "lastDispY"
  chosenDisp <- chosenCells[chosenCells[, "cellType"] %in% c(4, 3), , drop = FALSE]
  
  # Test
  if(testON == TRUE) { # GB_change P(sim)$testON
    expect_equivalent(NROW(chosenMat) + NROW(chosenDisp), length(unique(chosenCells[, "id"])))
  }
  
  if(NROW(chosenDisp) != 0) {
    dispersingIndDisp <- turtle(turtles = dispersingInd,
                                who = dispersingID[chosenDisp[, "id"]])
    dispersingInd <- NLset(turtles = dispersingInd, agents = dispersingIndDisp,
                           var = c("lastDispX", "lastDispY", "nMat"),
                           val = cbind(lastDispX = chosenDisp[, "pxcor"],
                                       lastDispY = chosenDisp[, "pycor"],
                                       nMat = 0))
  }
  
  # Movement
  chosenCellsCoords <- cbind(pxcor = chosenCells[, "pxcor"], pycor = chosenCells[, "pycor"])
  dispersingInd <- face(turtles = dispersingInd, agents2 = chosenCellsCoords)
  # Individuals that went back to their last known dispersal cell
  stop()
  if(NLcount(dispersingIndNMatMax) != 0){
    dispersingIndNMatMaxHead <- of(agents = dispersingIndNMatMax, var = "heading")
    headChoice <- sapply(dispersingIndNMatMaxHead,
                         function(x) {
                           which.min(abs(c(0, 45, 90, 135, 180, 225, 270, 315) - x))
                         })
    #dispersingInd <- NLset(turtles = dispersingInd, agents = dispersingIndNMatMax,
    #                       var = "heading", val = c(0, 45, 90, 135, 180, 225, 270, 315)[headChoice])
    dispersingInd@.Data[match(dispersingIndNMatMax@.Data[, "who"],
                              dispersingInd@.Data[, "who"]), "heading"] <- as.numeric(c(0, 45, 90, 135, 180, 225, 270, 315)[headChoice]) # faster
    
    # Reset dispersingIndNMatMax
    dispersingIndNMatMax <- noTurtles()
  }
  dispersingInd <- moveTo(turtles = dispersingInd, agents = chosenCellsCoords)
  # chosenCellsCoords are the dispersers new locations
  # Put + 1 on these cells because dispersers stepped on them
  sim$connectivityMap <- NLset(world = sim$connectivityMap, agents = chosenCellsCoords,
                               val = of(world = sim$connectivityMap, agents = chosenCellsCoords) + 1)
  
  # Test
  if(testON == TRUE) { # GB_change P(sim)$testON
    expect_equivalent(chosenCellsCoords, patchHere(world = sim$habitatMap,
                                                   turtles = dispersingInd))
    expect_equivalent(currCells[, c("pxcor", "pycor")], of(agents = dispersingInd,
                                                           var = c("prevX", "prevY")))
    expect_equivalent(patchHere(world = sim$habitatMap, turtles = dispersingInd)
                      [chosenCells[, "cellType"] %in% c(3, 4), ],
                      of(agents = dispersingInd, var = c("lastDispX", "lastDispY"))
                      [chosenCells[, "cellType"] %in% c(3, 4), ])
  }
  
  # Spatial mortality influenced by roads
  roadMort <- of(world = sim$roadMortMap, agents = chosenCellsCoords)
  # Add the correction factor for the dispersers
  deathRoad <- rbinom(n = length(roadMort), size = 1, prob = (roadMort / corrFactorDisp)) # GB_change P(sim)$corrFactorDisp
  # Do not kill the dispersers the first year of simulation
  # because all individuals from the initial population are dispersers
  ## INITIAL MORTALITY
  if(floor(time(sim))[1] == start(sim, "year")[1]){
    deathRoad <- rep(0, length(roadMort))
  }
  
  sim$nColl <- rbind(sim$nColl, data.frame(ncoll = sum(deathRoad), 
                                           time = floor(time(sim))[1]))
  deathRoad[roadMort == 1] <- 1 # mortality of 1 on the borders (roadMort == 1) needs to be forced
  deadWhoRoad <- dispersingID[deathRoad == 1]
  
  #stop()
  
  if(length(deadWhoRoad) != 0){
    sim$deadLynxColl[[time(sim, "year")[1]]] <- 
      turtleSet(sim$deadLynxColl[[time(sim, "year")[1]]], 
                turtle(turtles = dispersingInd, who = deadWhoRoad)) # add the new lynx dead by collisions
  }
  dispersingInd <- die(turtles = dispersingInd, who = deadWhoRoad)
  sim$aliveDispersingIndID <- dispersingInd@.Data[, "who"]
  sim$deadDisp[sim$deadDisp$time == floor(time(sim))[1], "nDispDeadColl"] <- 
    sim$deadDisp[sim$deadDisp$time == floor(time(sim))[1], "nDispDeadColl"] + length(deadWhoRoad)
  
  # Territory search
  disperser <- turtleSet(dispersingInd, nonDispersingInd)
  disperser <- sortOn(agents = disperser, var = "who")
  disperserID <- disperser@.Data[, "who"]
  sim$lynx <- turtleSet(disperser, nonDisperser)
  sim$lynx <- sortOn(agents = sim$lynx, var = "who")
  #stop()
  #source("C:/Users/gbal/Desktop/lynx.ibm/appendix_lynxIBM/c2.searchterritory.defunc.r")
  source("C:/Users/gbal/Desktop/lynx.ibm/appendix_lynxIBM/c1.searchterritory.func.r")
  start.searchTerr <- Sys.time()
  sim <- searchTerritory(sim)
  print(paste('searchTerr', start.searchTerr - Sys.time()))
  disperser <- turtle(turtles = sim$lynx, who = disperserID)
  nonDisperser <- turtle(turtles = sim$lynx, who = nonDisperserID)
  
} # end of if(NLcount(dispersingInd) != 0){

} # end of if(NLcount(disperser) != 0)
} # end of number of steps during the day
}



