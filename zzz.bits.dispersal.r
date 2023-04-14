# check of some headings change from Sarah's code
lynx.gb <- sim$lynx[,]
lynx.gb$heading
lynx.gb$heading %>% table
dispersingIndNMatMax.gb <- data.frame(who = 1150, heading = 90, prevX = 200, prevY = 200,
                                      breed = 'turtle', color = '#060702')
sapply(dispersingIndNMatMax.gb$heading,
       function(x) {
         which.min(abs(c(0, 45, 90, 135, 180, 225, 270, 315) - x))
       })

# return orderes alive
n.lynx.dummy <- 10
lynx.dummy.index <- sample.int(n.lynx.dummy, n.lynx.dummy, replace = FALSE)
lynx.dummy <- data.frame(
  ind = 1:n.lynx.dummy,
  who = letters[1:n.lynx.dummy],
  death = rbinom(n.lynx.dummy, 1, .3))[lynx.dummy.index, ]
who.order <- order(lynx.dummy$who)
who.alive.order <- who.order[lynx.dummy$death[who.order] == 0]
lynx.dummy[who.alive.order, ]
