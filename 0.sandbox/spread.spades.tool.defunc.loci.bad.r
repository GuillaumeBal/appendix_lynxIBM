if (getRversion() >= "3.1.0") {
  utils::globalVariables(c(".", ".I", "dists", "dup", "id", "indices", "initialLocus"))
}

require(namespace)
try(ns <- namespace::makeNamespace("spreadsDTInNamespace"))

landscape = 1; loci = NA_real_; spreadProb = 0.23; persistence = 0;
mask = NA; maxSize = 1e8L; directions = 8L; iterations = 1e6L;
lowMemory = NULL; # getOption("spades.lowMemory");
returnIndices = FALSE;
returnDistances = FALSE; mapID = NULL; id = FALSE; plot.it = FALSE;
spreadProbLater = NA_real_; spreadState = NA;
circle = FALSE; circleMaxRadius = NA_real_;
stopRule = NA; stopRuleBehavior = "includeRing"; allowOverlap = FALSE;
asymmetry = NA_real_; asymmetryAngle = NA_real_; quick = FALSE;
neighProbs = NULL; exactSizes = FALSE; relativeSpreadProb = FALSE


count.dim <- 1
dim.allcells.1 <- 1
while(dim.allcells.1 <= 96){
  x.picked <- sample.int(sim$habitatMap@maxPxcor, 1)
  y.picked <- sample.int(sim$habitatMap@maxPycor, 1)
  
  
  landscape = availCellsUpdatedRas
  loci = cellFromPxcorPycor(world = sim$habitatMap,
                            pxcor = x.picked,
                            pycor = y.picked)
  spreadProb = availCellsUpdatedRas 
  maxSize = terrSize
  returnIndices = TRUE
  quick = TRUE
  
  samInt <- dqrng::dqsample.int
  
  spreadStateExists <- is(spreadState, "data.table")
  spreadProbLaterExists <- TRUE
  
  ### should sanity check map extents
  
  loci <- as.integer(loci)
  
  msEqZero <- maxSize < 1
  initialLoci <- loci
  lenInitialLoci <- length(initialLoci)
  sequenceInitialLoci <- seq(lenInitialLoci)
  
  ncells <- as.integer(ncell(landscape))
  
  #browser(expr = exists("aaaaa"))
  allowOverlapOrReturnDistances <- allowOverlap | returnDistances
  useMatrixVersionSpreads <- allowOverlapOrReturnDistances | spreadStateExists
  
  # The experimental new spread function has some changes for speed. 1) The
  # bottleneck amazingly, was the creation of a new empty vector of length
  # ncell(landscape) ... it took >50% of the time of the spread function
  # when called 100,000s of times on a variety of spreadProb situations. 2) I
  # found that the only way to stop instantiating this was to have a
  # data.table object that uses reference semantics. 3) Put a simple, 1 column
  # data.table object into the SpaDES.tools namespace. It will contain the
  # former spreads object which was 0 everywhere the events hadn't spread
  # to, and a non-zero integer otherwise. 4) The function has to make sure that
  # it is "correct" on leaving the function. Two different cases: A) it
  # exits improperly --> action is delete this object; B) it exits correctly
  # --> action is to change all the values that were non-zero back to zero,
  # rather than delete the object. The whole point is to keep the object
  # intact after it has exited spread, so that it is available again
  # immediately for reuse.
  needEmptySpreads <- TRUE
  stNamespace <- asNamespace("SpaDES.tools")
  if (exists("spreadsDTInNamespace", envir = stNamespace)) {
    spreadsDT <- getFromNamespace("spreadsDTInNamespace", "SpaDES.tools")
    # set(spreadsDT, NULL, "spreads", 0L)
    # spreads <- spreadsDT$spreads
    if (identical(NROW(spreadsDT), ncells)) {
      needEmptySpreads <- FALSE
    }
  }
  if (needEmptySpreads) {
    spreadsDT <- data.table(spreads = vector("integer", ncells))
    set(spreadsDT, NULL, "spreads", 0L)
    # put the empty data.table into the SpaDES.tools namespace
    #assignInMyNamespace("spreadsDTInNamespace", spreadsDT)
    on.exit({assignInMyNamespace("spreadsDTInNamespace", integer())})
  }
  
  n <- 1L
  
  set(spreadsDT, loci, "spreads", 1L:length(loci))
  spreadsIndices <- unname(loci)
  #browser(expr = exists("aaaaa"))
  length(spreadsIndices) <- length(loci) * 100
  prevSpreadIndicesActiveLen <- length(loci)
  prevSpreadIndicesFullLen <- length(spreadsIndices)
  
  
  # Recycling maxSize as needed
  maxSize <- floor(maxSize)
  maxSize <- rep_len(maxSize, length(loci))
  size <- rep_len(1L, length(loci))
  
  #browser(expr = exists("aaaaa"))
  noMaxSize <- all(maxSize >= ncells) # will be used to omit testing for maxSize
  numNeighs <- NULL
  
  toColumn <- c("to", "indices")
  
  #browser(expr = exists("aaaaa"))
  # while there are active cells
  while (length(loci) & (n <= iterations)) {
    
    # identify neighbours
    potentials <- adj(landscape, loci, directions, pairs = TRUE)
    #stop()
    
    # keep only neighbours that have not been spread to yet
    # Keep only the ones where it hasn't been spread to yet
    keep <- spreadsDT$spreads[potentials[, 2L]] == 0L
    # keep <- spreads[potentials[, 2L]] == 0L
    potentials <- potentials[keep, , drop = FALSE]
    
    if (n == 2) {
      spreadProb <- spreadProbLater
    }
    
    # extract spreadProb values from spreadProb argument
    if(is.numeric(spreadProb)) {
      if (n == 1 & spreadProbLaterExists) {
        # need cell specific values
        spreadProbs <- rep(spreadProb, NROW(potentials))
        spreadProb <- spreadProbLater
      } else {
        if (length(spreadProb) > 1) {
          spreadProbs <- spreadProb[potentials[, 2L]]
        } else {
          spreadProbs <- rep(spreadProb, NROW(potentials))
        }
      }
    } else {
      # here for raster spreadProb
      if (n == 1 & spreadProbLaterExists) {
        # need cell specific values
        spreadProbs <- spreadProb[][potentials[, 2L]]
        spreadProb <- spreadProbLater
      } else {
        spreadProbs <- spreadProb[][potentials[, 2L]]
      }
    }
    
    if (anyNA(spreadProbs)) spreadProbs[is.na(spreadProbs)] <- 0
    
    randomSuccesses <- runifC(NROW(potentials)) <= spreadProbs
    potentials <- potentials[randomSuccesses, , drop = FALSE]
    
    # random ordering so not always same:
    lenPot <- NROW(potentials)
    
    # here is where allowOverlap and returnDistances are different ##### NOW OBSOLETE, I BELIEVE ELIOT March 2020
    potentials <- potentials[!duplicated(potentials[, 2L]), , drop = FALSE]
    
    # increment iteration
    n <- n + 1L
    
    # potentials can become zero because all active cells are edge cells
    #stop()
    if (length(potentials) > 0) {
      
      # If potentials has distances in it, it will be a numeric matrix; events should be integer
      events <- if (!is.integer(potentials)) as.integer(potentials[, 2L]) else potentials[, 2L]
      
      if (!noMaxSize) {
        #stop()
        len <- tabulate(spreadsDT$spreads[potentials[, 1L]], length(maxSize))
        size <- pmin(size + len, maxSize) ## Quick? and dirty. fast but loose (too flexible)
      }
      
      if (length(events) > 0) {
        if (id | returnIndices > 0 | relativeSpreadProb) {
          # give new cells, the id of the source cell
          set(spreadsDT, events, "spreads", spreadsDT$spreads[potentials[, 1L]])
        } #else {
        #set(spreadsDT, events, "spreads", n)
        #}
        curEventsLen <- length(events)
        addedIndices <- prevSpreadIndicesActiveLen + 1:curEventsLen
        
        # if (sum(curEventsLen, prevSpreadIndicesActiveLen) > prevSpreadIndicesFullLen) {
        #   length(spreadsIndices) <- (prevSpreadIndicesActiveLen + curEventsLen) * 2
        #   prevSpreadIndicesFullLen <- length(spreadsIndices)
        # }
        spreadsIndices[addedIndices] <- events
        prevSpreadIndicesActiveLen <- prevSpreadIndicesActiveLen + curEventsLen
        # spreadsIndices <- c(spreadsIndices, events)
        
      }
      
      # remove the cells from "events" that push it over maxSize
      #  There are some edge cases that aren't captured above ... identify them...
      if (length(maxSize) > 1L) {
        if (exists("whichID", inherits = FALSE)) {
          # must update toKeepSR in case that is a second reason to stop event
          if (exists("toKeepSR", inherits = FALSE)) {
            if (allowOverlapOrReturnDistances) {
              maxSizeKeep <- !(spreads[spreads[, "active"] == 1, "id"] %in% whichID)
              spreads <- spreads[c(rep(TRUE, sum(spreads[, "active"] == 0)), maxSizeKeep), ]
            } else {
              if (useMatrixVersionSpreads) {
                maxSizeKeep <- !spreads[events] %in% whichID
              } else {
                maxSizeKeep <- !spreadsDT$spreads[events] %in% whichID
              }
            }
            events <- events[maxSizeKeep]
            toKeepSR <- toKeepSR[maxSizeKeep]
          }
          rm(whichID)
        }
      } else {
        if (all(size >= maxSize)) {
          potentials <- potentials[0L,] # remove any potential cells, as size is met
          events <- NULL
        }
      }
      
      # Remove cells that were stopped by stopRule
      if (is.function(stopRule)) {
        if (exists("toKeepSR", inherits = FALSE)) {
          events <- events[toKeepSR]
          if (allowOverlapOrReturnDistances) {
            spreads[c(rep(TRUE, sum(spreads[, "active"] == 0)), !toKeepSR), "active"] <- 0
          }
          rm(toKeepSR)
        }
      }
    } else {
      # there are no potentials -- possibly from failed runif, or spreadProbs all 0
      events <- NULL
    }
    
    if (exactSizes) {
      if (all(get("numRetries", inherits = FALSE, envir = .pkgEnv) < 10)) {
        if (spreadStateExists) {
          tooSmall <- tabulate(spreads[, "id"], length(maxSize)) < maxSize
          inactive <- tabulate(spreads[spreads[, "active"] == 1, "id"], length(maxSize)) == 0
        } else {
          if (useMatrixVersionSpreads) {
            tooSmall <- tabulate(spreads, length(maxSize)) < maxSize
            inactive <- tabulate(spreads[events], length(maxSize)) == 0
          } else {
            tooSmall <- tabulate(spreadsDT$spreads, length(maxSize)) < maxSize
            inactive <- tabulate(spreadsDT$spreads[events], length(maxSize)) == 0
          }
        }
        
        # these are ones that are stuck ... i.e., too small, and inactive
        needPersist <- tooSmall & inactive
        needPersistJump <- TRUE
        if (any(needPersist)) {
          assign("numRetries", envir = .pkgEnv,
                 get("numRetries", inherits = FALSE, envir = .pkgEnv) + needPersist)
          
          if (spreadStateExists) {
            whSmallInactive <- which(tooSmall & inactive)
            spreadsSmallInactive <- spreads[spreads[, "id"] %in% whSmallInactive, , drop = FALSE]
            if (needPersistJump) {
              message("Jumping to new active location, up to 1000 m away")
              mmm <- rings(landscape, loci = spreadsSmallInactive[, "indices"],
                           maxRadius = 1000, minRadius = 1, returnIndices = TRUE)
              wh <- mmm[, list(whKeepLoci = resample(.I, 1)), by = id]$whKeepLoci
            } else {
              for (whSI in whSmallInactive) {
                wh <- which(spreads[, "id"] == whSI)
                wh <- tail(wh, 2) # pick last two ones from all inactive cells
                keepLoci <- spreads[wh, "indices"]
                events <- c(keepLoci, events)
                spreads[wh, "active"] <- 1
              }
            }
          } else {
            keepLoci <- spreads[loci] %in% which(tooSmall & inactive)
            events <- c(loci[keepLoci], events)
          }
        }
      }
    }
    
    # drop or keep loci
    
    loci <- NULL
    
    # new loci list for next while loop, concat of persistent and new events
    loci <- c(loci, events)
  } # end of while loop
  
  # Reset the base R seed so it is deterministic
  if (requireNamespace("dqrng", quietly = TRUE))
    set.seed(dqrng::dqsample.int(1e9, 1) + sample.int(1e9, 1))
  
  
  spreadsIndices <- spreadsIndices[1:prevSpreadIndicesActiveLen]
  
  
  # Convert the data back to raster
  
  wh <- 
    spreadsIndices
  
  wh <- wh[!(wh %in% potentials[,2L])]
  completed <- data.table(indices = wh, id = spreadsDT$spreads[wh], active = FALSE)
  
  active <- data.table(indices = as.integer(potentials[, 2L]),
                       id = spreadsDT$spreads[potentials[, 1L]],
                       active = TRUE)
  
  active <- data.table(indices = integer(0), id = integer(0), active = logical(0))
  
  if (returnIndices == 1) {
    #browser(expr = exists("aaaaa"))
    allCells <- rbindlist(list(completed, active)) # active first; next line will keep active
    initEventID <- allCells[indices %in% initialLoci, id]
    
    attr(initialLoci, ".match.hash") <- NULL # something in data.table put this
    dtToJoin <- data.table(id = sort(initEventID), initialLocus = initialLoci)
    
    setkeyv(dtToJoin, "id")
    setkeyv(allCells, "id")
    
    # tack on initialLoci
    allCells <- dtToJoin[allCells]
    
    #browser(expr = exists("aaaaa"))
    set(spreadsDT, allCells$indices, "spreads", 0L)
    # remove the previous on.exit which had the effect of deleting the contents
    #   completely on a failed `spread`. Here, we want to delete the previous
    #   on.exit --> allowing the object to stay intact, but with only zeros.
    on.exit()
    # spreadsIndices <- spreadsIndices[1:prevSpreadIndicesActiveLen]
    
  }  
  
  allCells
  dim.allcells.1 <- dim(allCells)[1]
  print(paste(count.dim, dim.allcells.1, sep = ' / '))
  count.dim <- count.dim + 1
} 
  
  
  