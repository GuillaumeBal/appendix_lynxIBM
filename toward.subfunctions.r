#outputs from towards
#315   0  45 270  90 225 180 135 135

radians.mov <- 
atan2(agents2[, 1] - agents[, 1], agents2[, 2] - agents[, 2])

RadToDegrees <- function(radian){
  degrees <- radian * (180 / pi)
  for(i in 1:length(degrees)){
    if(degrees[i]< 0){
      degrees[i] <- degrees[i] + 360
    }
  }
  return(degrees)
}

RadToDegrees(radians.mov)
