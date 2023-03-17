towards <- function(agents, agents2, world, torus = FALSE) {
  standardGeneric("towards")
})

#' @export
#' @rdname towards
setMethod(
  "towards",
  signature = c(agents = "matrix", agents2 = "matrix"),
  definition = function(agents, agents2, world, torus) {
    
    if (!inherits(agents, "agentMatrix") & !inherits(agents2, "agentMatrix")) {
      # patches to patches
      if (torus == FALSE) {
        heading <- deg(atan2(agents2[, 1] - agents[, 1], agents2[, 2] - agents[, 2]))
        # angles between -180 and 180
        heading[heading < 0] <- heading[heading < 0] + 360
      } else {
        if (missing(world)) {
          stop("A world must be provided as torus = TRUE")
        }
        
        if (NROW(agents2) == 1 & NROW(agents) != 1) {
          agents2 <- c(rep(agents2[, 1], NROW(agents)), rep(agents2[, 2], NROW(agents)))
          dim(agents2) <- c(NROW(agents), 2L)
        }
        if (NROW(agents) == 1 & NROW(agents2) != 1) {
          agents <- c(rep(agents[, 1], NROW(agents2)), rep(agents[, 2], NROW(agents2)))
          dim(agents) <- c(NROW(agents2), 2L)
        }
        
        # Need to create coordinates for "agents2" in a wrapped world
        # For all the 8 possibilities of wrapping (to the left, right, top, bottom and 4 corners)
        # Find the smallest distances across or around the world
        
        to1 <- cbind(agents2[, 1] - (world@extent@xmax - world@extent@xmin),
                     agents2[, 2] + (world@extent@ymax - world@extent@ymin))
        to2 <- cbind(agents2[, 1], agents2[, 2] + (world@extent@ymax - world@extent@ymin))
        to3 <- cbind(agents2[, 1] + (world@extent@xmax - world@extent@xmin),
                     agents2[, 2] + (world@extent@ymax - world@extent@ymin))
        to4 <- cbind(agents2[, 1] - (world@extent@xmax - world@extent@xmin), agents2[, 2])
        to5 <- cbind(agents2[, 1] + (world@extent@xmax - world@extent@xmin), agents2[, 2])
        to6 <- cbind(agents2[, 1] - (world@extent@xmax - world@extent@xmin),
                     agents2[, 2] - (world@extent@ymax - world@extent@ymin))
        to7 <- cbind(agents2[, 1], agents2[, 2] - (world@extent@ymax - world@extent@ymin))
        to8 <- cbind(agents2[, 1] + (world@extent@xmax - world@extent@xmin),
                     agents2[, 2] - (world@extent@ymax - world@extent@ymin))
        
        # All distances in a wrapped world
        distAgents2 <- pointDistance(p1 = agents, p2 = agents2, lonlat = FALSE, allpairs = FALSE)
        distTo1 <- pointDistance(p1 = agents, p2 = to1, lonlat = FALSE, allpairs = FALSE)
        distTo2 <- pointDistance(p1 = agents, p2 = to2, lonlat = FALSE, allpairs = FALSE)
        distTo3 <- pointDistance(p1 = agents, p2 = to3, lonlat = FALSE, allpairs = FALSE)
        distTo4 <- pointDistance(p1 = agents, p2 = to4, lonlat = FALSE, allpairs = FALSE)
        distTo5 <- pointDistance(p1 = agents, p2 = to5, lonlat = FALSE, allpairs = FALSE)
        distTo6 <- pointDistance(p1 = agents, p2 = to6, lonlat = FALSE, allpairs = FALSE)
        distTo7 <- pointDistance(p1 = agents, p2 = to7, lonlat = FALSE, allpairs = FALSE)
        distTo8 <- pointDistance(p1 = agents, p2 = to8, lonlat = FALSE, allpairs = FALSE)
        
        # Which distance is the minimum
        allDist <- cbind(distAgents2, distTo1, distTo2, distTo3, distTo4, distTo5,
                         distTo6, distTo7, distTo8)
        distMin <- apply(allDist, 1, min)
        
        toShortest <- agents2
        for (i in 1:NROW(agents)) {
          # All the possibilities for each agents (i.e., agents2 and the wrapped agents2)
          allToCoords <- rbind(agents2[i, ], to1[i, ], to2[i, ], to3[i, ], to4[i, ], to5[i, ],
                               to6[i, ], to7[i, ], to8[i, ])
          toShortest[i, ] <- allToCoords[match(distMin[i], allDist[i, ]), ]
          # if ties, take the first match (good because favor the non wrapped distances)
        }
        
        heading <- deg(atan2(toShortest[, 1] - agents[, 1], toShortest[, 2] - agents[, 2]))
        # angles between -180 and 180
        heading[heading < 0] <- heading[heading < 0] + 360
      }
      
    } else if (inherits(agents, "agentMatrix") & !inherits(agents2, "agentMatrix")) {
      # turtles to patches
      tCoords <- agents@.Data[, c("xcor", "ycor"), drop = FALSE]
      heading <- towards(agents = tCoords, agents2 = agents2, world = world, torus = torus)
      sameLoc <- tCoords[, 1] == agents2[, 1] & tCoords[, 2] == agents2[, 2]
      if (NROW(tCoords) == 1) {
        heading[sameLoc] <- agents@.Data[, "heading"]
      } else {
        heading[sameLoc] <- agents@.Data[, "heading"][sameLoc]
      }
      
    } else if (!inherits(agents, "agentMatrix") & inherits(agents2, "agentMatrix")) {
      # patches to turtles
      heading <- towards(agents = agents,
                         agents2 = agents2@.Data[, c("xcor", "ycor"), drop = FALSE],
                         world = world, torus = torus)
      
    } else if (inherits(agents, "agentMatrix") & inherits(agents2, "agentMatrix")) {
      # turtles to turtles
      t1Coords <- agents@.Data[, c("xcor", "ycor"), drop = FALSE]
      t2Coords <- agents2@.Data[, c("xcor", "ycor"), drop = FALSE]
      heading <- towards(agents = t1Coords, agents2 = t2Coords, world = world, torus = torus)
      sameLoc <- t1Coords[, 1] == t2Coords[, 1] & t1Coords[, 2] == t2Coords[, 2]
      if (NROW(t1Coords) == 1) {
        heading[sameLoc] <- agents@.Data[, "heading"]
      } else {
        heading[sameLoc] <- agents@.Data[, "heading"][sameLoc]
      }
    }
    
    return(heading)
  }
)