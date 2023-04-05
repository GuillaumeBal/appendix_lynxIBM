disp <- turtle(turtles = sim$lynx, who = sim$aliveDispersingIndID)

if(NLcount(disp) != 0) {
  
  # Female dispersers
  dispFem <- NLwith(agents = disp, var = "sex", val = "F")
  if(NLcount(dispFem) != 0) {
    dispFemID <- dispFem@.Data[, "who"]
    # Shuffle dispFemID so that it's not always the smallest IDs (i.e., the oldest individuals) that go first
    dispFemID <- dispFemID[sample(length(dispFemID))]
    
     searchingFemID <- dispFemID[1]
    #for(searchingFemID in dispFemID) {
      # To build territory = empty cells of breeding type (= 4)=
      searchingFemCell <- patchHere(world = sim$habitatMap, turtles = turtle(turtles = sim$lynx, who = searchingFemID))
      #searchingFemCellType <- of(world = sim$habitatMap, agents = searchingFemCell)
      searchingFemCellType <- sim$habitatMap[searchingFemCell[, 1], searchingFemCell[, 2]] # faster
      #searchingFemCellAvail <- of(world = sim$terrMap, agents = searchingFemCell)
      searchingFemCellAvail <- sim$terrMap[searchingFemCell[, 1], searchingFemCell[, 2]] # faster
      
      if(searchingFemCellType == 4 & is.na(searchingFemCellAvail)){ # current position
        #terrValues <- of(world = sim$terrMap, agents = patches(sim$terrMap))
        terrValues <- as.numeric(t(sim$terrMap@.Data)) # faster
        availCellsUpdatedRas <- sim$availCellsRas
        occupiedCells <- which(!is.na(terrValues))
        availCellsUpdatedRas[occupiedCells] <- 0
        
        # Spread from the female position (loci) to available contiguous cells
        # Find on which territory the female is
        terrSizeName <- of(world = sim$popDist, agents = patchHere(world = sim$popDist, turtles = turtle(turtles = sim$lynx, who = searchingFemID)))
        # And assign the territory size according to her position
        if(terrSizeName == 1){
          terrSize <- round(coreTerrSizeFAlps)
        }
        if(terrSizeName == 2){
          terrSize <- round(coreTerrSizeFJura)
        }
        if(terrSizeName == 3){
          terrSize <- round(coreTerrSizeFVosgesPalatinate)
        }
        if(terrSizeName == 4){
          terrSize <- round(coreTerrSizeFBlackForest)
        }
        terrDT <- spread(landscape = availCellsUpdatedRas,
                         loci = cellFromPxcorPycor(world = sim$habitatMap,
                                                   pxcor = searchingFemCell[, 1],
                                                   pycor = searchingFemCell[, 2]),
                         spreadProb = availCellsUpdatedRas, maxSize = terrSize, returnIndices = TRUE,
                         quick = TRUE)
        terrCells <- unique(terrDT$indices) # cells of the built territory
        if(length(terrCells) == terrSize) {
          newTerrCells <- PxcorPycorFromCell(world = sim$habitatMap, cellNum = terrCells)
          
          # Test
          if(PtestON == TRUE) { 
            expect_true(all(is.na(of(world = sim$terrMap, agents = newTerrCells))))
          }
          
          # Claim the territory
          sim$terrMap <- NLset(world = sim$terrMap, agents = newTerrCells, val = searchingFemID)
          # Mortality probability in the territory
          #probMortRdTerr <- of(world = sim$roadMortMap, agents = newTerrCells)
          probMortRdTerr <- sim$roadMortMap[newTerrCells[, 1], newTerrCells[, 2]] # faster
          sim$lynx <- NLset(turtles = sim$lynx, agents = turtle(turtles = sim$lynx, who = searchingFemID),
                            var = c("status", "rdMortTerr"),
                            val = cbind(status = "res", rdMortTerr = mean(probMortRdTerr)))
          # Save the data about the new residents
          if(length(searchingFemID) != 0){
            sim$resLynx[[time(sim, "year")[1]]] <- turtleSet(sim$resLynx[[time(sim, "year")[1]]], turtle(turtles = sim$lynx, who = searchingFemID))
            sim$timeRes <- rbind(sim$timeRes, data.frame(who = searchingFemID, year =  time(sim, "year")[1], time = sim$day))
          }
          
          # Male around to claim the female?
          # Check for a male on an area equal to the home range size (95 % kernel density) of the males in the population
          # Extract the distance to which look for a male as the radius of its home range size
          if(terrSizeName == 1){
            maxDistMale <- sqrt(terrSizeMAlps/pi)
          }
          if(terrSizeName == 2){
            maxDistMale <- sqrt(terrSizeMJura/pi)
          }
          if(terrSizeName == 3){
            maxDistMale <- sqrt(terrSizeMVosgesPalatinate/pi)
          }
          if(terrSizeName == 4){
            maxDistMale <- sqrt(terrSizeMBlackForest/pi)
          }
          neighbTerrCells <- NetLogoR::inRadius(agents = turtle(turtles = sim$lynx, who = searchingFemID), radius = maxDistMale,
                                                agents2 = patches(sim$habitatMap), world = sim$habitatMap, torus = FALSE)
          
          #neighbTerrCells <- unique(of(world = sim$terrMap, agents = neighbTerrCells))
          neighbTerrCells <- unique(sim$terrMap[neighbTerrCells[, 1], neighbTerrCells[, 2]]) # faster
          otherFemTerr <- neighbTerrCells[!is.na(neighbTerrCells) & neighbTerrCells != searchingFemID]
          if(length(otherFemTerr) != 0) {
            otherFem <- turtle(turtles = sim$lynx, who = otherFemTerr)
            infoOtherFem <- otherFem@.Data[, "maleID"]
            infoOtherFem <- infoOtherFem[!is.na(infoOtherFem)]
            otherMal <- turtle(turtles = sim$lynx, who = infoOtherFem)
            infoOtherMal <- otherMal@.Data[, "nFem"]
            
            if(length(infoOtherFem[infoOtherMal < 3]) != 0){
              # Calculate the distances between the new female resident and all the available males
              distFemaleMales <- NetLogoR::NLdist(agents = turtle(turtles = sim$lynx, who = searchingFemID),
                                                  agents2 = turtle(turtles = sim$lynx, who = infoOtherFem[infoOtherMal < 3]),
                                                  torus = FALSE)
              selectedMal <- ifelse(length(infoOtherFem[infoOtherMal < 3]) == 1,
                                    infoOtherFem[infoOtherMal < 3],
                                    # Select the closest male if there are several available
                                    infoOtherFem[infoOtherMal < 3][which.min(distFemaleMales)])
              
              # Associate the male to the female
              sim$lynx <- NLset(turtles = sim$lynx, agents = turtle(turtles = sim$lynx, who = searchingFemID),
                                var = "maleID", val = selectedMal)
              selectedMalInd <- turtle(turtles = sim$lynx, who = selectedMal)
              sim$lynx <- NLset(turtles = sim$lynx, agents = selectedMalInd, var = "nFem",
                                val = selectedMalInd@.Data[, "nFem"] + 1)
            }
          }
          
        } # end if(length(terrCells) == terrSize)
      } # end if(searchingFemCellType == 4 & is.na(searchingFemCellAvail))
    } # end for(searchingFemID in dispFemID)
  } # end if(NLcount(dispFem) != 0)
  
  # Male dispersers
  dispMal <- NLwith(agents = disp, var = "sex", val = "M")
  if(NLcount(dispMal) != 0) {
    dispMalID <- dispMal@.Data[, "who"]
    
    # Need to find a female territory
    cellMale <- patchHere(world = sim$habitatMap, turtles = dispMal)
    #cellMaleTerr <- of(world = sim$terrMap, agents = cellMale)
    cellMaleTerr <- sim$terrMap[cellMale[, 1], cellMale[, 2]] # faster
    dispMalIDwFem <- dispMalID[!is.na(cellMaleTerr)]
    
    # Shuffle dispMalIDwFem so that it's not always the smallest IDs (i.e., the oldest individuals) that go first
    dispMalIDwFem <- dispMalIDwFem[sample(length(dispMalIDwFem))]
    for(searchingMaleID in dispMalIDwFem) {
      
      whoFemEncountered <- cellMaleTerr[!is.na(cellMaleTerr)][
        dispMalID[!is.na(cellMaleTerr)] == searchingMaleID]
      femEncountered <- turtle(turtles = sim$lynx, who = whoFemEncountered)
      femEncounteredAvail <- femEncountered@.Data[, "maleID"]
      
      if(is.na(femEncounteredAvail)) {
        
        # Claim the female
        sim$lynx <- NLset(turtles = sim$lynx, agents = femEncountered, var = "maleID", val = searchingMaleID)
        maleClaiming <- turtle(turtles = sim$lynx, who = searchingMaleID)
        sim$lynx <- NLset(turtles = sim$lynx, agents = maleClaiming, var = c("status", "nFem"),
                          val = cbind(status =  "res", nFem = 1))
        # Save the data about the new residents
        if(length(maleClaiming) != 0){
          # There can be duplicates here as males becoming resident once can become disperser again if their female(s) die(s)
          # and so when they become resident again afterwards (in the same year) they are duplicated here.
          # However, turtleSet() only keep the first instance so there are no duplicates in the output resLynx, 
          # only the first time will be recorded.
          sim$resLynx[[time(sim, "year")[1]]] <- turtleSet(sim$resLynx[[time(sim, "year")[1]]], turtle(turtles = sim$lynx, who = of(agents = maleClaiming, var = "who")))
          sim$timeRes <- rbind(sim$timeRes, data.frame(who = of(agents = maleClaiming, var = "who"), year =  time(sim, "year")[1], time = sim$day))
        }
        
        # Females around to claim?
        # Check for females on an area equal to the male home range size (95 % kernel density) in the population
        # Extract the distance to which look for a male as the radius of its home range size
        terrSizeName <- of(world = sim$popDist, agents = patchHere(world = sim$popDist, turtles = turtle(turtles = sim$lynx, who = searchingMaleID)))
        if(terrSizeName == 1){
          maxDistMale <- sqrt(terrSizeMAlps/pi)
        }
        if(terrSizeName == 2){
          maxDistMale <- sqrt(terrSizeMJura/pi)
        }
        if(terrSizeName == 3){
          maxDistMale <- sqrt(terrSizeMVosgesPalatinate/pi)
        }
        if(terrSizeName == 4){
          maxDistMale <- sqrt(terrSizeMBlackForest/pi)
        }
        neighbTerrCells <- NetLogoR::inRadius(agents = turtle(turtles = sim$lynx, who = searchingMaleID), radius = maxDistMale,
                                              agents2 = patches(sim$habitatMap), world = sim$habitatMap, torus = FALSE)
        
        neighbTerrCells <- unique(sim$terrMap[neighbTerrCells[, 1], neighbTerrCells[, 2]]) # faster
        otherFemTerr <- neighbTerrCells[!is.na(neighbTerrCells) & neighbTerrCells != whoFemEncountered]
        if(length(otherFemTerr) != 0) {
          otherFem <- turtle(turtles = sim$lynx, who = otherFemTerr)
          infoOtherFem <- otherFem@.Data[, "maleID"]
          if(length(otherFemTerr[is.na(infoOtherFem)]) != 0){
            # Calculate the distances between the new male resident and all the available females 
            distMaleFemales <- NetLogoR::NLdist(agents = turtle(turtles = sim$lynx, who = searchingMaleID),
                                                agents2 = turtle(turtles = sim$lynx, who = otherFemTerr[is.na(infoOtherFem)]),
                                                torus = FALSE)
            
            # Select up to 2 females without male
            selectedFem <- ifelse(length(otherFemTerr[is.na(infoOtherFem)]) >= 2,
                                  # Select the closest females if there are more than 2
                                  otherFemTerr[is.na(infoOtherFem)][which.minn(distMaleFemales, n = 2)],
                                  otherFemTerr[is.na(infoOtherFem)])
            
            
            # Claim the female(s)
            sim$lynx <- NLset(turtles = sim$lynx, agents = turtle(turtles = sim$lynx, who = selectedFem),
                              var = "maleID", val = searchingMaleID)
            sim$lynx <- NLset(turtles = sim$lynx, agents = maleClaiming, var = "nFem",
                              val = 1 + length(selectedFem))
          }
        }
        
      } # end if(is.na(femEncounteredAvail))
    } # end for(searchingMaleID in dispMalIDwFem)
  } # end if(NLcount(dispMal) != 0)
  
  # Test
  if(testON == TRUE) {
    terrNumber <- of(world = sim$terrMap, agents = patches(sim$terrMap))
    terrNumber <- terrNumber[!is.na(terrNumber)]
    expect_true(all(table(terrNumber) >= min(c(coreTerrSizeFAlps, coreTerrSizeFJura, coreTerrSizeFVosgesPalatinate,
                                               coreTerrSizeFBlackForest))))
    expect_true(all(table(terrNumber) <= max(c(coreTerrSizeFAlps, coreTerrSizeFJura, coreTerrSizeFVosgesPalatinate, 
                                               coreTerrSizeFBlackForest))))
    infoPop <- sim$lynx@.Data[, c("who", "maleID", "nFem"), drop = FALSE]
    nFemPerMal <- table(infoPop[, "maleID"])
    expect_equivalent(as.numeric(nFemPerMal),
                      infoPop[infoPop[, "who"] %in% as.numeric(names(nFemPerMal)), "nFem"])
    expect_true(all(sim$lynx@.Data[, "nFem"] >= 0 & sim$lynx@.Data[, "nFem"] <= 3))
    terrNumTerrMap <- unique(of(world = sim$terrMap, agents = patches(sim$terrMap)))
    terrNumLynx <- of(agents = NLwith(agents = NLwith(agents = sim$lynx, var = "sex", val = "F"),
                                      var = "status", val = "res"), var = "who")
    expect_true(all(terrNumTerrMap[!is.na(terrNumTerrMap)] %in% terrNumLynx))
  }
}


