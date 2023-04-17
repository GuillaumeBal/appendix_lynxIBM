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


landscape = availCellsUpdatedRas
loci = cellFromPxcorPycor(world = sim$habitatMap,
                          pxcor = 520,#searchingFemCell[, 1],
                          pycor = 10)#searchingFemCell[, 2])
spreadProb = availCellsUpdatedRas 
maxSize = terrSize
returnIndices = TRUE
quick = TRUE

if (!is.null(neighProbs)) {
  if (isTRUE(allowOverlap))
    stop("Can't use neighProbs and allowOverlap = TRUE together")
}
if (requireNamespace("dqrng", quietly = TRUE)) {
  samInt <- dqrng::dqsample.int
  # set dqrng seed from base state
  dqrng::dqset.seed(sample.int(1e9, 2))
} else {
  samInt <- sample.int
}

if (!is.null(mapID)) {
  warning("mapID is deprecated, use id")
  id <- mapID
}
if (!quick) {
  allowedRules <- c("includePixel", "excludePixel", "includeRing", "excludeRing")
  if (!any(stopRuleBehavior %in% allowedRules))
    stop("stopRuleBehaviour must be one of \"",
         paste(allowedRules, collapse = "\", \""), "\".")
}
if (isTRUE(lowMemory)) {
  stop("lowMemory is no longer supported due to removal of ffbase from CRAN.")
  # requireNamespace("ff", quietly = TRUE)
  # requireNamespace("ffbase", quietly = TRUE)
}

spreadStateExists <- is(spreadState, "data.table")
spreadProbLaterExists <- TRUE

if (!is(spreadProbLater, "Raster")) {
  if (anyNA(spreadProbLater)) {
    spreadProbLaterExists <- FALSE
    spreadProbLater <- spreadProb
  }
}

### should sanity check map extents
if (any(is.na(loci)))  {
  # start it in the centre cell, if there is no spreadState
  if (!spreadStateExists)
    loci <- middlePixel(landscape) #(nrow(landscape) / 2L + 0.5) * ncol(landscape)
}
if (!is.integer(loci)) {
  loci <- as.integer(loci)
}
if (!quick) {
  dupLoci <- duplicated(loci)
  if (any(duplicated(loci))) {
    message("duplicate initial loci are provided")
    # loci <- loci[dupLoci]
  }
}

if (length(loci) == 0) stop("No loci. Nothing to do")

if (any(!is.na(maxSize))) {
  msEqZero <- maxSize < 1
  if (any(msEqZero)) {
    loci <- loci[!msEqZero]
    maxSize <- maxSize[!msEqZero]
  }
}

if (spreadStateExists) {
  keepers <- spreadState$active == TRUE
  loci <- initialActiveCells <- spreadState[keepers, indices]
  initialLoci <- unique(spreadState$initialLocus)
} else {
  initialLoci <- loci
}
lenInitialLoci <- length(initialLoci)
sequenceInitialLoci <- seq(lenInitialLoci)

# Check for probabilities
if (!quick) {
  if (is(spreadProbLater, "RasterLayer") | is(spreadProb, "Rasterlayer")) {
    if ((minValue(spreadProb) > 1L) || (maxValue(spreadProb) < 0L) ||
        (maxValue(spreadProb) > 1L) || (minValue(spreadProb) < 0L)) {
      relativeSpreadProb <- TRUE
    }
    if (spreadProbLaterExists)
      if (((minValue(spreadProbLater) > 1L) || (maxValue(spreadProbLater) < 0L) ||
           (maxValue(spreadProbLater) > 1L) || (minValue(spreadProbLater) < 0L))) {
        relativeSpreadProb <- TRUE
      }
  } else {
    if (!all(inRange(na.omit(spreadProb)))) {
      relativeSpreadProb <- TRUE
      stop("spreadProb is not a probability")
    }
    if (spreadProbLaterExists) {
      relativeSpreadProb <- TRUE
      if (!all(inRange(na.omit(spreadProbLater)))) stop("spreadProbLater is not a probability")
    }
  }
}

ncells <- as.integer(ncell(landscape))

#browser(expr = exists("aaaaa"))
allowOverlapOrReturnDistances <- allowOverlap | returnDistances
useMatrixVersionSpreads <- allowOverlapOrReturnDistances | spreadStateExists
if (useMatrixVersionSpreads) {
  if (spreadStateExists) {
    spreads <- as.matrix(spreadState[, list(initialLocus, indices, id, active)])
  } else {
    spreads <- cbind(initialLocus = initialLoci, indices = initialLoci,
                     id = 1:length(loci), active = 1)
  }
} else {
  if (!is.null(lowMemory)) {
    message("lowMemory argument is now deprecated; using standard spread.")
  }
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
}

n <- 1L

# circle needs directions to be 8
if (circle | !is.na(asymmetry)) {
  if (circle) directions <- 8L # only required for circle
  initialLociXY <- cbind(id = seq_along(initialLoci), xyFromCell(landscape, initialLoci))
  id <- TRUE
  if (allowOverlapOrReturnDistances) {
    spreads <- cbind(spreads, dists = 0)
  }
}

# determine ... variables
# otherVars <- list(...)
# anyList <- unlist(lapply(otherVars, is.list))
# 
# if (any(anyList)) {
#   otherVarsLists <- unlist(unname(otherVars), recursive = FALSE)
#   otherVars[anyList] <- NULL
#   otherVars <- append(otherVars, otherVarsLists)
# }

# check validity of stopRule
if (is.function(stopRule)) {
  id <- TRUE
  stopRuleObjs <- names(formals(stopRule))
  if (!quick) {
    if (any(is.na(match(stopRuleObjs,
                        c("id", "landscape", "cells", names(otherVars)))))) {
      stop("Arguments in stopRule not valid.\n",
           "The function definition must be a function of built-in options,",
           " (id, landscape, or cells) or user supplied variables.",
           " If user supplied, the variables",
           " must be passed as named vectors, or lists or data.frames.",
           " See examples.")
    }
  }
  landRasNeeded <- any(stopRuleObjs == "landscape")
  colNamesPotentials <- c("id", "landscape"[landRasNeeded], "cells", "prev")
  argNames <- c(colNamesPotentials, names(otherVars))
  whArgs <- match(names(formals(stopRule)), argNames)
  
  # Raster indexing is slow. If there is are Rasters submitted with the stopRule
  #  then this will convert them to vectors first. Clearly, this will have
  #  memory consequences if the Rasters are on disk, but spread is optimized for speed
  rasters <- unlist(lapply(otherVars[names(otherVars)], function(x) is(x, "Raster")))
  if (any(rasters)) {
    for (i in 1:which(rasters)) {
      otherVars[[names(rasters[i])]] <- otherVars[[names(rasters[i])]][]
    }
  }
  landRas <- landscape[] # For speed
}

if (!allowOverlap && !returnDistances) {
  if (id | returnIndices > 0 | relativeSpreadProb) {
    if (!spreadStateExists) {
      set(spreadsDT, loci, "spreads", 1L:length(loci))
      ##DT spreads[loci] <- 1L:length(loci)
      # give values to spreads vector at initialLoci
    }
  } else {
    spreadsDT$spreads[loci] <- n
  }
  spreadsIndices <- unname(loci)
  #browser(expr = exists("aaaaa"))
  length(spreadsIndices) <- length(loci) * 100
  prevSpreadIndicesActiveLen <- length(loci)
  prevSpreadIndicesFullLen <- length(spreadsIndices)
}

# Convert mask and NAs to 0 on the spreadProb Raster
if (is(spreadProb, "Raster")) {
  # convert NA to 0s
  #isNASpreadProb <- is.na(spreadProb[])
  # if (anyNA(spreadProb[])) {
  #   isNASpreadProb <- is.na(spreadProb[])
  #   spreadProb[isNASpreadProb] <- 0L
  # }
} else if (is.numeric(spreadProb)) {
  # Translate numeric spreadProb into a Raster, if there is a mask
  if (is(mask, "Raster")) {
    spreadProb <- raster(extent(landscape), res = res(landscape), vals = spreadProb)
  }
}

# Convert mask and NAs to 0 on the spreadProbLater Raster
if (is(spreadProbLater, "Raster")) {
} else if (is.numeric(spreadProbLater)) {
  # Translate numeric spreadProbLater into a Raster, if there is a mask
  if (is(mask, "Raster")) {
    spreadProbLater <- raster(extent(landscape), res = res(landscape), vals = spreadProbLater)
  }
}

# Mask spreadProbLater and spreadProb
if (is(mask, "Raster")) {
  spreadProbLater[mask[] == 1L] <- 0L
  spreadProb[mask[] == 1L] <- 0L
}

if (spreadStateExists) {
  if (allowOverlapOrReturnDistances) {
    stop("Using spreadState with either allowOverlap = TRUE",
         " or returnDistances = TRUE is not implemented")
  } else {
    if (sum(colnames(spreadState) %in% c("indices", "id", "active", "initialLocus")) != 4) {
      stop("spreadState must have at least columns: ",
           "indices, id, active, and initialLocus.")
    }
  }
}

if (!quick) {
  if (any(loci > ncells)) stop("loci indices are not on landscape")
}

# Recycling maxSize as needed
if (any(!is.na(maxSize))) {
  if (!is.integer(maxSize)) maxSize <- floor(maxSize)
  if (spreadStateExists) {
    sizeAll <- spreadState[, list(len = .N), by = id]
    size <- c(sizeAll[, len])
  } else {
    maxSize <- rep_len(maxSize, length(loci))
    size <- rep_len(1L, length(loci))
  }
} else {
  maxSize <- ncells
  size <- length(loci)
}

#browser(expr = exists("aaaaa"))
noMaxSize <- all(maxSize >= ncells) # will be used to omit testing for maxSize
if (is.null(neighProbs)) {
  numNeighs <- NULL
}

# if (!exists("numRetries", envir = .pkgEnv)) {
#   assign("numRetries", rep(0, lenInitialLoci), envir = .pkgEnv)
# }

toColumn <- c("to", "indices")

#browser(expr = exists("aaaaa"))
# while there are active cells
while (length(loci) & (n <= iterations)) {
  
  if (!is.null(neighProbs)) {
    numNeighs <- if (is.list(neighProbs)) {
      unlist(lapply(neighProbs, function(x) {
        sample.int(length(x), size = 1, replace = TRUE, prob = x)
      }))
    } else {
      sample.int(length(neighProbs), size = length(loci), replace = TRUE, prob = neighProbs)
    }
  }
  
  # identify neighbours
  if (useMatrixVersionSpreads) {
    whActive <- spreads[, "active"] == 1 # spreads carries over
    potentials <- adj(landscape, loci, directions, pairs = TRUE,
                      id = spreads[whActive, "id"])#, numNeighs = numNeighs)
    spreads[whActive, "active"] <- 0
    potentials <- cbind(potentials, active = 1)
  } else {
    if (id | returnIndices > 0 | circle | relativeSpreadProb | !is.null(neighProbs)) {
      potentials <- adj(landscape, loci, directions, pairs = TRUE)
    } else {
      # must pad the first column of potentials
      newAdj <- adj(landscape, loci, directions, pairs = FALSE)
      potentials <- cbind(NA_integer_, newAdj)
    }
  }
  
  if (circle) {
    potentials <- cbind(potentials, dists = 0)
  }
  
  # keep only neighbours that have not been spread to yet
  if (useMatrixVersionSpreads) {
    # data.table version is faster for potentials > 2000 or so
    if (NROW(potentials) > 2000) {
      spreadsDT <- as.data.table(spreads)
      potentialsDT <- as.data.table(potentials)
      potentialsDT[, initialLocus := initialLoci[potentialsDT$id]]
      colnamesPot <- colnames(potentialsDT)
      whIL <- which(colnamesPot == "initialLocus")
      whFrom <- which(colnamesPot == "from")
      setcolorder(potentialsDT,
                  c(colnamesPot[whIL], colnamesPot[-c(whIL, whFrom)], colnamesPot[whFrom]))
      setnames(potentialsDT, old = "to", new = "indices")
      newPot <- potentialsDT[!spreadsDT, on = c("id", "indices")]
      potentials <- as.matrix(newPot)
    } else {
      potentials <- cbind(initialLocus = initialLoci[potentials[, "id"]], potentials)
      colnames(potentials)[which(colnames(potentials) == "to")] <- "indices"
      colnamesPot <- colnames(potentials)
      whIL <- which(colnamesPot == "initialLocus")
      whFrom <- which(colnamesPot == "from")
      potentials <- potentials[,c(colnamesPot[whIL],
                                  colnamesPot[-c(whIL, whFrom)],
                                  colnamesPot[whFrom])]
      
      # These next lines including the lapply are the rate limiting
      #   step and it has been heavily worked to speed it up March 31, 2020
      seq2 <- sequenceInitialLoci[sequenceInitialLoci %in% potentials[,"id"]]
      out <- lapply(seq2, function(ind) {
        hasID <- potentials[, "id"] == ind
        po <- potentials[hasID, ]
        hasID2 <- spreads[, "id"] == ind
        inds <- spreads[hasID2, "indices"]
        vals <- po[, 2L] %in% inds
        po[!vals,]
      })
      potentials <- do.call(rbind, out)
    }
  } else {
    # Keep only the ones where it hasn't been spread to yet
    keep <- spreadsDT$spreads[potentials[, 2L]] == 0L
    # keep <- spreads[potentials[, 2L]] == 0L
    potentials <- potentials[keep, , drop = FALSE]
  }
  
  if (n == 2) {
    spreadProb <- spreadProbLater
  }
  
  # extract spreadProb values from spreadProb argument
  if (is.numeric(spreadProb)) {
    if (!(length(spreadProb) == 1 || length(spreadProb) == ncell(landscape)))
      stop("spreadProb must be length 1 or length ncell(landscape), or a raster")
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
  
  if (!is.na(asymmetry)) {
    if (allowOverlapOrReturnDistances) {
      a <- cbind(id = potentials[, 3L], to = potentials[, 2L],
                 xyFromCell(landscape, potentials[, 2L]))
    } else {
      if (useMatrixVersionSpreads) {
        a <- cbind(id = spreads[potentials[, 1L]], to = potentials[, 2L],
                   xyFromCell(landscape, potentials[, 2L]))
      } else {
        a <- cbind(id = spreadsDT$spreads[potentials[, 1L]], to = potentials[, 2L],
                   xyFromCell(landscape, potentials[, 2L]))
      }
    }
    d <- directionFromEachPoint(from = initialLociXY, to = a)
    newSpreadProbExtremes <- (spreadProb[] * 2) / (asymmetry + 1) * c(1, asymmetry)
    angleQuality <- (cos(d[, "angles"] - rad(asymmetryAngle)) + 1) / 2
    spreadProbs <- newSpreadProbExtremes[1] + (angleQuality * diff(newSpreadProbExtremes))
    spreadProbs <- spreadProbs - diff(c(spreadProb[], mean(spreadProbs)))
  }
  
  if (!is.null(neighProbs) | relativeSpreadProb) {
    aaa <- split(seq_along(potentials[, toColumn[spreadStateExists + 1]]),
                 potentials[, "from"]);
    if (length(aaa) != length(numNeighs)) {
      activeCellContinue <- loci %in% unique(potentials[, "from"])
      numNeighs <- numNeighs[activeCellContinue]
    }
    
    tmpA <- unlist(lapply(aaa, length))
    tmpB <- which(tmpA < numNeighs)
    if (length(tmpB) > 0)
      numNeighs[tmpB] <- unname(tmpA[tmpB])
    
    if (relativeSpreadProb) {
      rescaledProbs <- tapply(spreadProbs, potentials[, "from"], function(x) {
        x / sum(x, na.rm = TRUE)
      }, simplify = FALSE)
      neighIndexToKeep <- unlist(lapply(seq_along(aaa), function(x)
        resample(aaa[[x]], size = numNeighs[x], prob = rescaledProbs[[x]])))
    } else {
      neighIndexToKeep <- unlist(lapply(seq_along(aaa), function(x)
        resample(aaa[[x]], size = numNeighs[x])))
    }
    potentials <- potentials[neighIndexToKeep, , drop = FALSE]
    spreadProbs <- spreadProbs[neighIndexToKeep]
    spreadProbs[spreadProbs > 0] <- 1
  }
  
  randomSuccesses <- runifC(NROW(potentials)) <= spreadProbs
  potentials <- potentials[randomSuccesses, , drop = FALSE]
  
  # random ordering so not always same:
  lenPot <- NROW(potentials)
  if (lenPot) {
    reorderVals <- samInt(lenPot)
    potentials <- potentials[reorderVals, , drop = FALSE]
  }
  if (!allowOverlap) {
    # here is where allowOverlap and returnDistances are different ##### NOW OBSOLETE, I BELIEVE ELIOT March 2020
    potentials <- potentials[!duplicated(potentials[, 2L]), , drop = FALSE]
  } else {
    pots <- potentials[, c("id", "indices"), drop = FALSE]
    potentials <- potentials[!duplicated(pots), , drop = FALSE]
  }
  
  # increment iteration
  n <- n + 1L
  
  # potentials can become zero because all active cells are edge cells
  if (length(potentials) > 0) {
    # implement circle
    if (!missing(circle)) {
      if (circle) {
        if (allowOverlapOrReturnDistances) {
          a <- cbind(potentials, xyFromCell(landscape, potentials[, 2L]))
        } else {
          #browser(expr = exists("aaaaa"))
          if (useMatrixVersionSpreads) {
            a <- cbind(potentials, id = spreads[potentials[, "from"]],
                       xyFromCell(landscape, potentials[, "to"]))
          } else {
            a <- cbind(potentials, id = spreadsDT$spreads[potentials[, "from"]],
                       xyFromCell(landscape, potentials[, "to"]))
          }
        }
        # need to remove dists column because distanceFromEachPoint, adds one back
        a <- a[, !(colnames(a) %in% c("dists")), drop = FALSE]
        # need 3 columns, id, x, y in both initialLociXY and a
        d <- distanceFromEachPoint(initialLociXY, a, angles = asymmetry) # d is sorted
        cMR <- (n - 1) * res(landscape)[1]
        if (!any(is.na(circleMaxRadius))) {
          # don't bother proceeding if circleMaxRadius is larger than current iteration
          if (any(circleMaxRadius <= ((n - 1) * res(landscape)[1]))) {
            if (length(circleMaxRadius) > 1) {
              # if it is a vector of values
              cMR <- circleMaxRadius[d[, "id"]]
            } else {
              cMR <- circleMaxRadius
            }
          }
        }
        potentials <- d[, !(colnames(d) %in% c("x", "y")), drop = FALSE]
        potentials <- potentials[(d[, "dists"] %<=% cMR), , drop = FALSE]
      }
    }
    
    # If potentials has distances in it, it will be a numeric matrix; events should be integer
    events <- if (!is.integer(potentials)) as.integer(potentials[, 2L]) else potentials[, 2L]
    
    if (!noMaxSize) {
      if (useMatrixVersionSpreads) {
        len <- tabulate(potentials[, 3L], length(maxSize))
      } else {
        # actually interested in potential[,2L], but they don't have values yet..
        #  can use their source
        len <- tabulate(spreadsDT$spreads[potentials[, 1L]], length(maxSize))
      }
      if (any((size + len) > maxSize & size <= maxSize)) {
        whichID <- which(size + len > maxSize)
        
        # remove some active cells, if more than maxSize
        toRm <- (size + len)[whichID] - maxSize[whichID]
        
        for (i in 1:length(whichID)) {
          if (useMatrixVersionSpreads) {
            thisID <- which(potentials[, 3L] == whichID[i])
          } else {
            thisID <- which(spreadsDT$spreads[potentials[, 1L]] == whichID[i])
          }
          
          # some unusual cases where there are none on the spreads. Unsure how this occurs
          if (length(thisID))
            potentials <- potentials[-resample(thisID, toRm[i]), , drop = FALSE]
        }
        events <- as.integer(potentials[, 2L])
      }
      size <- pmin(size + len, maxSize) ## Quick? and dirty. fast but loose (too flexible)
    }
    
    # Implement stopRule section
    if (is.function(stopRule) & length(events) > 0) {
      if (allowOverlapOrReturnDistances) {
        prevCells <- cbind(
          id = spreads[, "id"],
          landscape = if (landRasNeeded) landRas[spreads[, "indices"]] else NULL,
          cells = spreads[, "indices"], prev = 1)
        eventCells <- cbind(id = potentials[, "id"],
                            landscape = if (landRasNeeded) landRas[events] else NULL,
                            cells = events, prev = 0)
      } else {
        whgtZero <- spreadsIndices
        #browser(expr = exists("aaaaa"))
        if (useMatrixVersionSpreads) {
          prevCells <- cbind(id = spreads[whgtZero],
                             landscape = if (landRasNeeded) landRas[whgtZero] else NULL,
                             cells = whgtZero, prev = 1)
          eventCells <- cbind(id = spreads[potentials[, 1L]],
                              landscape = if (landRasNeeded) landRas[potentials[, 2L]] else NULL,
                              cells = potentials[, 2L], prev = 0)
        } else {
          prevCells <- cbind(id = spreadsDT$spreads[whgtZero],
                             landscape = if (landRasNeeded) landRas[whgtZero] else NULL,
                             cells = whgtZero, prev = 1)
          eventCells <- cbind(id = spreadsDT$spreads[potentials[, 1L]],
                              landscape = if (landRasNeeded) landRas[potentials[, 2L]] else NULL,
                              cells = potentials[, 2L], prev = 0)
        }
      }
      if (circle) {
        prevCells <- cbind(prevCells, dist = NA)
        eventCells <- cbind(eventCells, dist = potentials[, "dists"])
      }
      # don't need to continue doing ids that are not active
      tmp <- rbind(prevCells[prevCells[, "id"] %in% unique(eventCells[, "id"]), ], eventCells)
      
      ids <- unique(tmp[, "id"])
      
      shouldStopList <- lapply(ids, function(id) {
        shortTmp <- tmp[tmp[, "id"] == id, ]
        args <- append(list(id = id),
                       lapply(colNamesPotentials[-1], function(j) shortTmp[, j]))
        names(args) <- colNamesPotentials
        args <- append(args, otherVars)
        do.call(stopRule, args[whArgs])
      })
      if (any(lapply(shouldStopList, length) > 1))
        stop("stopRule does not return a length-one logical.",
             " Perhaps stopRule need indexing by cells or id?")
      
      shouldStop <- unlist(shouldStopList)
      
      names(shouldStop) <- ids
      
      if (any(shouldStop)) {
        if (stopRuleBehavior != "includeRing") {
          if (stopRuleBehavior != "excludeRing") {
            whStop <- as.numeric(names(shouldStop)[shouldStop])
            whStopAll <- tmp[, "id"] %in% whStop
            tmp2 <- tmp[whStopAll, ]
            
            whStopEvents <- eventCells[, "id"] %in% whStop
            
            # If an event needs to stop, then must identify which cells are included
            out <- lapply(whStop, function(id) {
              tmp3 <- tmp2[tmp2[, "id"] == id, ]
              newOnes <- tmp3[, "prev"] == 0
              ord <- seq_along(newOnes)
              
              # because of undesired behaviour of sample when length(x) == 1
              if (sum(newOnes) > 1) {
                ord[newOnes] <- sample(ord[newOnes])
                if (circle) ord[newOnes] <- ord[newOnes][order(tmp3[ord[newOnes], "dist"])]
                tmp3 <- tmp3[ord, ]
              }
              startLen <- sum(!newOnes)
              addIncr <- 1
              done <- FALSE
              args <- append(list(id = id),
                             lapply(colNamesPotentials[-1], function(j) {
                               tmp3[1:startLen, j]
                             })) # instead of as.data.frame
              names(args) <- colNamesPotentials
              args <- append(args, otherVars)
              argsSeq <- seq_along(colNamesPotentials[-1]) + 1
              
              while (!done) {
                args[argsSeq] <- lapply(colNamesPotentials[-1], function(j) {
                  unname(c(args[[j]], tmp3[(startLen + addIncr), j]))
                }) # instead of as.data.frame
                done <- do.call(stopRule, args[whArgs])
                addIncr <- addIncr + 1
              }
              if (stopRuleBehavior == "excludePixel") addIncr <- addIncr - 1
              firstInd <- startLen + addIncr
              lastInd <- NROW(tmp3)
              sequ <- if (firstInd > lastInd) 0 else firstInd:lastInd
              tmp3[sequ, , drop = FALSE]
            })
            eventRm <- do.call(rbind, out)[, "cells"]
            cellsKeep <- !(potentials[, 2L] %in% eventRm)
          } else {
            cellsKeep <- rep(FALSE, NROW(potentials))
          }
          potentials <- potentials[cellsKeep, , drop = FALSE]
          events <- as.integer(potentials[, 2L])
          eventCells <- eventCells[cellsKeep, , drop = FALSE]
        }
        toKeepSR <- !(eventCells[, "id"] %in% as.numeric(names(which((shouldStop)))))
      }
    }
    
    if (length(events) > 0) {
      # place new value at new cells that became active
      if (useMatrixVersionSpreads) {
        fromCol <- colnames(potentials) == "from"
        
        spreads <- rbind(spreads, potentials[, !fromCol])
        if ((returnDistances | spreadStateExists) & !allowOverlap) {
          # 2nd place where allowOverlap and returnDistances differ
          notDups <- !duplicated(spreads[, "indices"])
          nrSpreads <- NROW(spreads)
          nrPotentials <- NROW(potentials)
          notDupsEvents <- notDups[-(1:(nrSpreads - nrPotentials))]
          spreads <- spreads[notDups, , drop = FALSE]
          events <- events[notDupsEvents]
        }
      } else {
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
        # spreadsIndices <- c(spreadsIndices, events)
      }
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
  if (is.na(persistence) | persistence == 0L) {
    loci <- NULL
  } else {
    if (inRange(persistence)) {
      loci <- loci[runif(length(loci)) <= persistence]
    } else {
      # here is were we would handle methods for raster* or functions
      stop("Unsupported type: persistence")
    }
  }
  
  if (plot.it) {
    if (n == 2 & !spreadStateExists) clearPlot()
    if (allowOverlapOrReturnDistances) {
      spreadsDT <- data.table(spreads);
      hab2 <- landscape;
      hab2[] <- 0;
      pixVal <- spreadsDT[, sum(id), by = indices]
      hab2[pixVal$indices] <- pixVal$V1;
      Plot(hab2, legendRange = c(0, sum(seq_along(initialLoci))))
    } else {
      plotCur <- raster(landscape)
      plotCur <- setValues(plotCur, spreads)
      Plot(plotCur)
    }
  }
  
  # new loci list for next while loop, concat of persistent and new events
  loci <- c(loci, events)
} # end of while loop

# Reset the base R seed so it is deterministic
if (requireNamespace("dqrng", quietly = TRUE))
  set.seed(dqrng::dqsample.int(1e9, 1) + sample.int(1e9, 1))

if (!allowOverlap & !returnDistances) {
  spreadsIndices <- spreadsIndices[1:prevSpreadIndicesActiveLen]
}

# Convert the data back to raster
if (!allowOverlap & !returnDistances & !spreadStateExists) {
  # if (lowMemory) {
  #   wh <- as.ram(ffwhich(spreads, spreads > 0))
  #   if (returnIndices > 0) {
  #     wh <- wh[!(wh %in% potentials[,2L])]
  #     completed <- data.table(indices = wh, id = spreads[wh], active = FALSE)
  #     if (NROW(potentials) > 0) {
  #       active <- data.table(indices = potentials[, 2L],
  #                            id = spreads[potentials[, 1L]],
  #                            active = TRUE)
  #     } else {
  #       active <- data.table(indices = numeric(0), id = numeric(0),
  #                            active = logical(0))
  #     }
  #   }
  # } else {
  wh <- if (spreadStateExists) {
    c(spreadState[!keepers]$indices, spreadsIndices)
  } else {
    spreadsIndices
  }
  if (returnIndices > 0) {
    # wh already contains the potentials for next iteration -- these should be not duplicated
    #   inside "completed"
    wh <- wh[!(wh %in% potentials[,2L])]
    completed <- data.table(indices = wh, id = spreadsDT$spreads[wh], active = FALSE)
    if (NROW(potentials) > 0) {
      active <- data.table(indices = as.integer(potentials[, 2L]),
                           id = spreadsDT$spreads[potentials[, 1L]],
                           active = TRUE)
    } else {
      active <- data.table(indices = integer(0), id = integer(0), active = logical(0))
    }
  }
}

if (returnIndices == 1) {
  if (useMatrixVersionSpreads) {
    keepCols <- c(3, 1, 2, 4)
    if (circle) keepCols <- c(keepCols, 5)
    
    # change column order to match non allowOverlap
    allCells <- as.data.table(spreads[, keepCols, drop = FALSE])
    set(allCells, NULL, j = "active", as.logical(allCells$active))
    # setkeyv(allCells, "id")
  } else {
    #browser(expr = exists("aaaaa"))
    allCells <- rbindlist(list(completed, active)) # active first; next line will keep active
    if (spreadStateExists) {
      initEventID <- unique(spreadState$id)
    } else {
      initEventID <- allCells[indices %in% initialLoci, id]
    }
    if (!all(is.na(initialLoci))) {
      attr(initialLoci, ".match.hash") <- NULL # something in data.table put this
      dtToJoin <- data.table(id = sort(initEventID), initialLocus = initialLoci)
    } else {
      dtToJoin <- data.table(id = numeric(0), initialLocus = numeric(0))
    }
    setkeyv(dtToJoin, "id")
    setkeyv(allCells, "id")
    
    # tack on initialLoci
    allCells <- dtToJoin[allCells]
  }
  allCells[]
  if (exactSizes)
    if (exists("numRetries", envir = .pkgEnv)) {
      if (sum(allCells$active) == 0) rm("numRetries", envir = .pkgEnv)
    }
  if  (!(useMatrixVersionSpreads)) {
    #browser(expr = exists("aaaaa"))
    set(spreadsDT, allCells$indices, "spreads", 0L)
    # remove the previous on.exit which had the effect of deleting the contents
    #   completely on a failed `spread`. Here, we want to delete the previous
    #   on.exit --> allowing the object to stay intact, but with only zeros.
    on.exit()
    # spreadsIndices <- spreadsIndices[1:prevSpreadIndicesActiveLen]
  }
  
}  

allCells
