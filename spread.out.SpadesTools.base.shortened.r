if (getRversion() >= "3.1.0") {
  utils::globalVariables(c(".", ".I", "dists", "dup", "id", "indices", "initialLocus"))
}

borrowedResample <- utils::getFromNamespace("resample", "SpaDES.tools")
#body(borrowedfun)

# resampleZeroProof <- function(spreadProbHas0, x, n, prob) {
#   if (spreadProbHas0) {
#     sm <- sum(prob, na.rm = TRUE)
#     if (sum(prob > 0) <= n) {
#       integer()
#     } else {
#       resample(x, n, prob = prob / sm)
#     }
#   } else resample(x, n, prob = prob / sum(prob, na.rm = TRUE))
# }


spread.short.gb <- 
  function(landscape, loci = NA_real_, spreadProb = 0.23, persistence = 0,
           mask = NA, maxSize = 1e8L, directions = 8L, iterations = 1e6L,
           lowMemory = NULL, # getOption("spades.lowMemory"),
           returnIndices = FALSE,
           returnDistances = FALSE, mapID = NULL, id = FALSE, plot.it = FALSE,
           spreadProbLater = NA_real_, spreadState = NA,
           circle = FALSE, circleMaxRadius = NA_real_,
           stopRule = NA, stopRuleBehavior = "includeRing", allowOverlap = FALSE,
           asymmetry = NA_real_, asymmetryAngle = NA_real_, quick = FALSE,
           neighProbs = NULL, exactSizes = FALSE, relativeSpreadProb = FALSE) {
    
    samInt <- sample.int
    spreadStateExists <- is(spreadState, "data.table")
    spreadProbLaterExists <- TRUE
    
    # essential if loop GB !!!!!!!!!!!!!!!!!!! 
    if (!is(spreadProbLater, "Raster")) {
      if (anyNA(spreadProbLater)) {
        spreadProbLaterExists <- FALSE
        spreadProbLater <- spreadProb
      }
    }
    
    loci <- as.integer(loci)
    initialLoci <- loci
    lenInitialLoci <- length(initialLoci)
    sequenceInitialLoci <- seq(lenInitialLoci)
    ncells <- as.integer(ncell(landscape))
    
    #browser(expr = exists("aaaaa"))
    allowOverlapOrReturnDistances <- allowOverlap | returnDistances
    useMatrixVersionSpreads <- allowOverlapOrReturnDistances | spreadStateExists
    
    needEmptySpreads <- TRUE
    
    spreadsDT <- data.table(spreads = vector("integer", ncells))
    set(spreadsDT, NULL, "spreads", 0L)
    # put the empty data.table into the SpaDES.tools namespace
    assign("spreadsDTInNamespace", spreadsDT)
    on.exit({assign("spreadsDTInNamespace", integer())})
    
    n <- 1L
    
    spreadsDT$spreads[loci] <- n
    
    spreadsIndices <- unname(loci)
    length(spreadsIndices) <- length(loci) * 100
    prevSpreadIndicesActiveLen <- length(loci)
    prevSpreadIndicesFullLen <- length(spreadsIndices)
    
    # Recycling maxSize as needed
    if (any(!is.na(maxSize))) {
      maxSize <- rep_len(maxSize, length(loci))
      size <- rep_len(1L, length(loci))
    } else {
      maxSize <- ncells
      size <- length(loci)
    }
    
    #browser(expr = exists("aaaaa"))
    noMaxSize <- all(maxSize >= ncells) # will be used to omit testing for maxSize
    if (is.null(neighProbs)) {
      numNeighs <- NULL
    }
    
    
    toColumn <- c("to", "indices")
    
    while (length(loci) & (n <= iterations)) {
      
      # identify neighbours
      potentials <- adj(landscape, loci, directions, pairs = TRUE)
      # Keep only the ones where it hasn't been spread to yet
      keep <- spreadsDT$spreads[potentials[, 2L]] == 0L
      # keep <- spreads[potentials[, 2L]] == 0L
      potentials <- potentials[keep, , drop = FALSE]
      
      spreadProbs <- spreadProb[potentials[, 2L]]
      if (anyNA(spreadProbs)) spreadProbs[is.na(spreadProbs)] <- 0
      
      randomSuccesses <- runifC(NROW(potentials)) <= spreadProbs
      potentials <- potentials[randomSuccesses, , drop = FALSE]
      
      # random ordering so not always same:
      lenPot <- NROW(potentials)
      
      reorderVals <- samInt(lenPot)
      potentials <- potentials[reorderVals, , drop = FALSE]
      
      
      # here is where allowOverlap and returnDistances are different ##### NOW OBSOLETE, I BELIEVE ELIOT March 2020
      potentials <- potentials[!duplicated(potentials[, 2L]), , drop = FALSE]
      
      # increment iteration
      n <- n + 1L
      
      # potentials can become zero because all active cells are edge cells
      if (length(potentials) > 0) {
        
        # If potentials has distances in it, it will be a numeric matrix; events should be integer
        events <- if (!is.integer(potentials)) as.integer(potentials[, 2L]) else potentials[, 2L]
        
        if (!noMaxSize) {
          #if (useMatrixVersionSpreads) {
          # len <- tabulate(potentials[, 3L], length(maxSize))
          #} else {
          # actually interested in potential[,2L], but they don't have values yet..
          #  can use their source
          len <- tabulate(spreadsDT$spreads[potentials[, 1L]], length(maxSize))
          #}
          if (any((size + len) > maxSize & size <= maxSize)) {
            whichID <- which(size + len > maxSize)
            
            # remove some active cells, if more than maxSize
            toRm <- (size + len)[whichID] - maxSize[whichID]
            
            for (i in 1:length(whichID)) {
              
              thisID <- which(spreadsDT$spreads[potentials[, 1L]] == whichID[i])
              
              # some unusual cases where there are none on the spreads. Unsure how this occurs
              if (length(thisID))
                potentials <- potentials[-borrowedResample(thisID, toRm[i]), , drop = FALSE]
            }
            events <- as.integer(potentials[, 2L])
          }
          size <- pmin(size + len, maxSize) ## Quick? and dirty. fast but loose (too flexible)
        }
        
        if (length(events) > 0) {
          if (id | returnIndices > 0 | relativeSpreadProb) {
            # give new cells, the id of the source cell
            set(spreadsDT, events, "spreads", spreadsDT$spreads[potentials[, 1L]])
          } else {
            set(spreadsDT, events, "spreads", n)
          }
          curEventsLen <- length(events)
          addedIndices <- prevSpreadIndicesActiveLen + 1:curEventsLen
          
          if (sum(curEventsLen, prevSpreadIndicesActiveLen) > prevSpreadIndicesFullLen) {
            length(spreadsIndices) <- (prevSpreadIndicesActiveLen + curEventsLen) * 2
            prevSpreadIndicesFullLen <- length(spreadsIndices)
          }
          spreadsIndices[addedIndices] <- events
          prevSpreadIndicesActiveLen <- prevSpreadIndicesActiveLen + curEventsLen
        }
        
        # remove the cells from "events" that push it over maxSize
        #  There are some edge cases that aren't captured above ... identify them...
        if (length(maxSize) > 1L) {
          potentials <- potentials[0L,] # remove any potential cells, as size is met
          events <- NULL
        }
        
      } else {
        # there are no potentials -- possibly from failed runif, or spreadProbs all 0
        events <- NULL
      }
      
      # drop or keep loci
      # if (is.na(persistence) | persistence == 0L) {
      loci <- NULL
      # new loci list for next while loop, concat of persistent and new events
      loci <- c(loci, events)
    } # end of while loop
    
    # Reset the base R seed so it is deterministic
    spreadsIndices <- spreadsIndices[1:prevSpreadIndicesActiveLen]
    
    # Convert the data back to raster
    if (!allowOverlap & !returnDistances & !spreadStateExists) {
      wh <- #if (spreadStateExists) {
        #   c(spreadState[!keepers]$indices, spreadsIndices)
        # } else {
        spreadsIndices
      # }
      if (returnIndices > 0) {
        # wh already contains the potentials for next iteration -- these should be not duplicated
        #   inside "completed"
        wh <- wh[!(wh %in% potentials[,2L])]
        completed <- data.table(indices = wh, id = spreadsDT$spreads[wh], active = FALSE)
        active <- data.table(indices = integer(0), id = integer(0), active = logical(0))
      }
    }
    
    if (returnIndices == 1) {
      
      allCells <- rbindlist(list(completed, active)) # active first; next line will keep active
      
      initEventID <- allCells[indices %in% initialLoci, id]
      
      attr(initialLoci, ".match.hash") <- NULL # something in data.table put this
      dtToJoin <- data.table(id = sort(initEventID), initialLocus = initialLoci)
      
      setkeyv(dtToJoin, "id")
      setkeyv(allCells, "id")
      
      # tack on initialLoci
      allCells <- dtToJoin[allCells]
      
      return(allCells)
    }
    
  }











