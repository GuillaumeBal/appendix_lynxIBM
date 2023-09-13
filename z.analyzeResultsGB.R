sex.dich <- c("M", "F")

movePopLongDT <- list()

for(s in sex.dich){
  
  movePop <- cbind(repSim = rep(1:nSim, each = lastYear), year = rep(1:lastYear, nSim), 
                   AtoA = rep(0, lastYear*nSim), AtoJ = rep(0, lastYear*nSim), AtoBF = rep(0, lastYear*nSim), AtoVP = rep(0, lastYear*nSim),
                   JtoJ = rep(0, lastYear*nSim), JtoA = rep(0, lastYear*nSim), JtoBF = rep(0, lastYear*nSim), JtoVP = rep(0, lastYear*nSim),
                   BFtoBF = rep(0, lastYear*nSim), BFtoA = rep(0, lastYear*nSim), BFtoJ = rep(0, lastYear*nSim), BFtoVP = rep(0, lastYear*nSim),
                   VPtoVP = rep(0, lastYear*nSim), VPtoA = rep(0, lastYear*nSim), VPtoJ = rep(0, lastYear*nSim), VPtoBF = rep(0, lastYear*nSim))
  
  for(i in 1:length(listSim)){ # for each simulation run
    load(paste0(pathFiles, "/", listSim[i]))
    
    for(y in 1:lastYear){
      
      # Identify the residents in the different populations
      resAlps <- turtlesOn(world = lynxIBMrun$popDist, turtles = lynxIBMrun$resLynx[[y]] %>% `[`(.$sex == s),
                           agents = NLwith(agents = patches(lynxIBMrun$popDist), world = lynxIBMrun$popDist, val = 1))
      resJura <- turtlesOn(world = lynxIBMrun$popDist, turtles = lynxIBMrun$resLynx[[y]] %>% `[`(.$sex == s),
                           agents = NLwith(agents = patches(lynxIBMrun$popDist), world = lynxIBMrun$popDist, val = 2))
      resVP <- turtlesOn(world = lynxIBMrun$popDist, turtles = lynxIBMrun$resLynx[[y]] %>% `[`(.$sex == s),
                         agents = NLwith(agents = patches(lynxIBMrun$popDist), world = lynxIBMrun$popDist, val = 3))
      resBF <- turtlesOn(world = lynxIBMrun$popDist, turtles = lynxIBMrun$resLynx[[y]] %>% `[`(.$sex == s),
                         agents = NLwith(agents = patches(lynxIBMrun$popDist), world = lynxIBMrun$popDist, val = 4))
      
      # Where these residents were from?
      if(NLcount(resAlps) != 0){
        whoResAlps <- of(agents = resAlps, var = "pop")
        movePop[movePop[, "year"] == y & movePop[, "repSim"] == i, "AtoA"] <- length(which(whoResAlps == "Alps")) # established in the Alps but coming from the Jura
        movePop[movePop[, "year"] == y & movePop[, "repSim"] == i, "JtoA"] <- length(which(whoResAlps == "Jura")) # established in the Alps but coming from the Jura
        movePop[movePop[, "year"] == y & movePop[, "repSim"] == i, "VPtoA"] <- length(which(whoResAlps == "Vosges-Palatinate")) # established in the Alps but coming from the Vosges-Palatinate
        movePop[movePop[, "year"] == y & movePop[, "repSim"] == i, "BFtoA"] <- length(which(whoResAlps == "BlackForest")) # established in the Alps but coming from the Black Forest
      }
      if(NLcount(resJura) != 0){
        whoResJura <- of(agents = resJura, var = "pop")
        movePop[movePop[, "year"] == y & movePop[, "repSim"] == i, "JtoJ"] <- length(which(whoResJura == "Jura")) 
        movePop[movePop[, "year"] == y & movePop[, "repSim"] == i, "AtoJ"] <- length(which(whoResJura == "Alps")) # established in the Jura but coming from the Alps
        movePop[movePop[, "year"] == y & movePop[, "repSim"] == i, "VPtoJ"] <- length(which(whoResJura == "Vosges-Palatinate")) # established in the Jura but coming from the Vosges-Palatinate
        movePop[movePop[, "year"] == y & movePop[, "repSim"] == i, "BFtoJ"] <- length(which(whoResJura == "BlackForest")) # established in the Jura but coming from the Black Forest
      }
      if(NLcount(resVP)){
        whoResVP <- of(agents = resVP, var = "pop")
        movePop[movePop[, "year"] == y & movePop[, "repSim"] == i, "VPtoVP"] <- length(which(whoResVP == "Vosges-Palatinate")) 
        movePop[movePop[, "year"] == y & movePop[, "repSim"] == i, "AtoVP"] <- length(which(whoResVP == "Alps")) # established in the Vosges-Palatinate but coming from the Alps
        movePop[movePop[, "year"] == y & movePop[, "repSim"] == i, "JtoVP"] <- length(which(whoResVP == "Jura")) # established in the Vosges-Palatinate but coming from the Jura
        movePop[movePop[, "year"] == y & movePop[, "repSim"] == i, "BFtoVP"] <- length(which(whoResVP == "BlackForest")) # established in the Vosges-Palatinate but coming from the Black Forest
      }
      if(NLcount(resBF)){
        whoResBF <- of(agents = resBF, var = "pop")
        movePop[movePop[, "year"] == y & movePop[, "repSim"] == i, "BFtoBF"] <- length(which(whoResBF == "BlackForest")) 
        movePop[movePop[, "year"] == y & movePop[, "repSim"] == i, "AtoBF"] <- length(which(whoResBF == "Alps")) # established in the Black Forest but coming from the Alps
        movePop[movePop[, "year"] == y & movePop[, "repSim"] == i, "JtoBF"] <- length(which(whoResBF == "Jura")) # established in the Black Forest but coming from the Jura
        movePop[movePop[, "year"] == y & movePop[, "repSim"] == i, "VPtoBF"] <- length(which(whoResBF == "Vosges-Palatinate")) # established in the Black Forest but coming from the Vosges-Palatinate
      }
    }
    print(i)
  }
  
  # Summarize the data
  movePopLongDT[[s]] <- melt(setDT(as.data.frame(movePop)), id = c("repSim", "year"))
  # Compute the cumulative sum of the movement over the years per simulation and per variable (=pair of populations with a direction)
  movePopLongDT[[s]][, Cum.Sum := cumsum(value), by=list(repSim, variable)] 
  
}

movePopLongDT[['F']] %<>% as.data.frame()
movePopLongDT[['F']][ , c('from', 'to')] <- movePopLongDT[['F']]$variable %>% as.character() %>%
  strsplit(x = ., split = 'to') %>% do.call(., what = 'rbind') 

movePopLongDT[['M']] %<>% as.data.frame()
movePopLongDT[['M']][ , c('from', 'to')] <- movePopLongDT[['M']]$variable %>% as.character() %>%
  strsplit(x = ., split = 'to') %>% do.call(., what = 'rbind') 

tapply(movePopLongDT[['F']]$value, INDEX = list(movePopLongDT[['F']]$from, movePopLongDT[['F']]$to), FUN = sum) %>% prop.table(., margin = 1) %>% round(., 3)
tapply(movePopLongDT[['M']]$value, INDEX = list(movePopLongDT[['M']]$from, movePopLongDT[['M']]$to), FUN = sum) %>% prop.table(., margin = 1) %>% round(., 3)

# survival juv region ================================================================




# stuff found while looking at a crs issue #####################################

proj_db <- system.file("proj/proj.db", package = "sf")
crs_table <- sf::read_sf(proj_db, "crs_view")

crs_table[crs_table$name %>% grep(., pattern = 'MERC'), ]
#habMapSpaDES@crs
habMapSpaDES
habMapSpaDES %>% class

