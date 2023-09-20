move.save[["F"]] %>% dim
move.save[["M"]] %>% dim

move.save[["M"]] %>% summary

check.corr <- (move.save[["F"]] %>% summary %>% colnames()) %in% 
  (move.save[["M"]] %>% summary %>% colnames())

female.names <- (move.save[["F"]] %>% summary %>% colnames())

female.names[!check.corr]

paste0(stats.to.fill %>% names, '.TO.', POP.CURR) %in% colnames(movePop)




n <- 170; 
team <- letters[factor(rep_len(1:3, n), levels = 1:5)]
year <- sample.int(n = 10, size = n, replace = TRUE)
victory <- rep(1, n)

tapply(victory, INDEX = list(team, year), FUN = sum)

table(fac)
tapply(fac, )
tapply(1:n, fac, sum)
tapply(1:n, fac, sum, default = 0) # maybe more desirable
tapply(1:n, fac, sum, simplify = FALSE)
tapply(1:n, fac, range)
tapply(1:n, fac, quantile)
tapply(1:n, fac, length) ## NA's
tapply(1:n, fac, length, default = 0) # == taby
