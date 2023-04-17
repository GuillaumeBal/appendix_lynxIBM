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
loci = as.integer(cellFromPxcorPycor(world = sim$habitatMap,
                                     pxcor = x.picked,
                                     pycor = y.picked)
                  )
spreadProb = availCellsUpdatedRas 
maxSize = terrSize
returnIndices = TRUE
quick = TRUE


samInt <- dqrng::dqsample.int
# set dqrng seed from base state
#dqrng::dqset.seed(sample.int(1e9, 2))

spreadStateExists <- is(spreadState, "data.table")
spreadProbLaterExists <- TRUE

### should sanity check map extents

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


# ALL THESE LINES MAKE A spreadsDT$spreads, A VECTOR WITH ONLY ZERO SIZE OF HABITAT TYPE MATRICES
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
maxSize <- ncells
size <- length(loci)


#browser(expr = exists("aaaaa"))
noMaxSize <- all(maxSize >= ncells) # will be used to omit testing for maxSize
numNeighs <- NULL


# if (!exists("numRetries", envir = .pkgEnv)) {
#   assign("numRetries", rep(0, lenInitialLoci), envir = .pkgEnv)
# }

toColumn <- c("to", "indices")

#browser(expr = exists("aaaaa"))
# while there are active cells
while (length(loci) & (n <= iterations)) {
  
  # identify neighbours
  potentials <- adj(landscape, loci, directions, pairs = TRUE)
  
  
  # keep only neighbours that have not been spread to yet
  
  # Keep only the ones where it hasn't been spread to yet
  keep <- spreadsDT$spreads[potentials[, 2L]] == 0L
  # keep <- spreads[potentials[, 2L]] == 0L
  potentials <- potentials[keep, , drop = FALSE]
  
  if (n == 2) {
    spreadProb <- spreadProbLater
  }
  
  # extract spreadProb values from spreadProb argument
  # here for raster spreadProb
  if (n == 1 & spreadProbLaterExists) {
    # need cell specific values
    spreadProbs <- spreadProb[][potentials[, 2L]]
    spreadProb <- spreadProbLater
  } else {
    spreadProbs <- spreadProb[][potentials[, 2L]]
  }
  
  #if (anyNA(spreadProbs)) spreadProbs[is.na(spreadProbs)] <- 0
  
  randomSuccesses <- runifC(NROW(potentials)) <= spreadProbs
  potentials <- potentials[randomSuccesses, , drop = FALSE]
  
  # random ordering so not always same:
  lenPot <- NROW(potentials)
  
  # here is where allowOverlap and returnDistances are different ##### NOW OBSOLETE, I BELIEVE ELIOT March 2020
  potentials <- potentials[!duplicated(potentials[, 2L]), , drop = FALSE]
  
  # increment iteration
  n <- n + 1L
  
  # potentials can become zero because all active cells are edge cells
  
  # there are no potentials -- possibly from failed runif, or spreadProbs all 0
  events <- NULL
  
  # drop or keep loci
  loci <- NULL
  
  # new loci list for next while loop, concat of persistent and new events
  loci <- c(loci, events)
} # end of while loop

# Reset the base R seed so it is deterministic
#set.seed(dqrng::dqsample.int(1e9, 1) + sample.int(1e9, 1))

spreadsIndices <- spreadsIndices[1:prevSpreadIndicesActiveLen]


# Convert the data back to raster

wh <- spreadsIndices
active <- data.table(indices = integer(0), id = integer(0), active = logical(0))

#browser(expr = exists("aaaaa"))
allCells <- rbindlist(list(completed, active)) # active first; next line will keep active

initEventID <- allCells[indices %in% initialLoci, id]

attr(initialLoci, ".match.hash") <- NULL # something in data.table put this
dtToJoin <- data.table(id = sort(initEventID), initialLocus = initialLoci)

setkeyv(dtToJoin, "id")
setkeyv(allCells, "id")

# tack on initialLoci
allCells <- dtToJoin[allCells]

allCells
