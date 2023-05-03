# map to look at dispersing

hab.map.with.lynx <- sim$habitatMap
hab.map.with.lynx@.Data[lynx.gb$lastDispY, lynx.gb$lastDispX] <- 50
hab.map.with.lynx@.Data[outputs.loop[[1]]$lynx_lastDispY, outputs.loop[[1]]$lynx_lastDispX] <- 200

plot(hab.map.with.lynx)
