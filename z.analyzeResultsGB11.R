sex.dich <- c("M", "F")
pop <- c(
  'auver.fr',  
  'forez.fr',
  'morvan.fr',
  'alpsS.fr',
  'jura.fr',
  'vosges.fr',
  'jura.ch',#7
  'plat.ch',
  'alps.ch',
  'blackf.de',
  'palat.de'
)
n.pop <- length(pop)
pop.comb <- 
  combn(pop, m = 2, FUN = NULL, simplify = FALSE) %>% do.call(., what = rbind) %>% 
  apply(., MARGIN = 1, FUN = paste, collapse = '.TO.', simplify = FALSE) %>% unlist

# a few things to create list to feed in loop
movePopLongDT <- list()
matrix.mov <- matrix(0, nrow = length(rep(0, lastYear*nSim)), ncol = length(pop.comb))
df.mov <- matrix.mov %>% as.data.frame() %>% `colnames<-`(pop.comb)

# renaming the regions based on ours
pop.dist.gb <- raster2world(raster::raster("map.11.pop.tif"))
pop.dist.gb %>% plot

for(s in sex.dich){
  #s <- "M"
  
  movePop <- cbind(repSim = rep(1:nSim, each = lastYear), year = rep(1:lastYear, nSim), 
                   df.mov)
  
  for(i in 1:length(listSim)){ # for each simulation run
    
    load(paste0(pathFiles, "/", listSim[i]))
    
    # some gn pop corrections
    lynxIBMrun$popDist <- pop.dist.gb
    
    lynxIBMrun.gb <- lynxIBMrun
    lynxIBMrun.gb$resLynx %>% sapply(class)
    lynxIBMrun.gb
    
    for(l in 1:lastYear){
      if(length(lynxIBMrun.gb$resLynx[[l]][,'pop']) != 0)
        lynxIBMrun.gb$resLynx[[l]][,'pop'] <- 
          pop[lynxIBMrun.gb$popDist[cbind(lynxIBMrun.gb$resLynx[[l]][,'lastDispY'],
                                          lynxIBMrun.gb$resLynx[[l]][,'lastDispX'])]] 
    }
    
    for(y in 1:lastYear){
      #y <- 1
      
      # Identify the residents in the different populations
      for(p in 1:n.pop){
        assign(paste0("RES", pop[p]),
               value =turtlesOn(world = lynxIBMrun.gb$popDist, turtles = lynxIBMrun.gb$resLynx[[y]] %>% `[`(.$sex == s),
                                agents = NLwith(agents = patches(lynxIBMrun.gb$popDist), world = lynxIBMrun.gb$popDist, val = p))
        )
      }
      
      # resAlps <- turtlesOn(world = lynxIBMrun$popDist, turtles = lynxIBMrun$resLynx[[y]] %>% `[`(.$sex == s),
      #                      agents = NLwith(agents = patches(lynxIBMrun$popDist), world = lynxIBMrun$popDist, val = 1))
      # resJura <- turtlesOn(world = lynxIBMrun$popDist, turtles = lynxIBMrun$resLynx[[y]] %>% `[`(.$sex == s),
      #                      agents = NLwith(agents = patches(lynxIBMrun$popDist), world = lynxIBMrun$popDist, val = 2))
      # resVP <- turtlesOn(world = lynxIBMrun$popDist, turtles = lynxIBMrun$resLynx[[y]] %>% `[`(.$sex == s),
      #                    agents = NLwith(agents = patches(lynxIBMrun$popDist), world = lynxIBMrun$popDist, val = 3))
      # resBF <- turtlesOn(world = lynxIBMrun$popDist, turtles = lynxIBMrun$resLynx[[y]] %>% `[`(.$sex == s),
      #                    agents = NLwith(agents = patches(lynxIBMrun$popDist), world = lynxIBMrun$popDist, val = 4))
      
      for(n in 1:n.pop){
        if(NLcount(get(paste0("RES", pop[p]))) != 0){
          POP.CURR <- pop[p]
          RES.CURR <- get(paste0("RES", POP.CURR))
          WHO.CURR <- of(agents = RES.CURR, var = "pop")
          stats.to.fill <- table(factor(whoResAlps, level = pop))
          movePop[movePop[, "year"] == y & movePop[, "repSim"] == i, paste0(pop, '.TO.', POP.CURR)] <- stats.to.fill
        }
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

