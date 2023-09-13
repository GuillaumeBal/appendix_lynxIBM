issue <- 41 + 1
issue.range <- seq(issue - 5, issue + 5, 1)
cell.index.pb.range <- cell.index[issue.range]
issue.range.XY <- CellNumtoRowCol(cell.index.pb.range - 1, sim$availCellsRas %>% as.matrix)

plot(sim$availCellsRas)
points(x = issue.range.XY$x_coords - 1, 
       y = issue.range.XY$y_coords - 1, 
       col = 'red', lwd = 2)



for(i in 1:length(cell.index.pb.range)){
  outputs.spread <- try(
    spreadGB(
      availCellsMat = sim$availCellsRas %>% as.matrix,
      XCoordinate = issue.range.XY$x_coords[i],#cell.map.max - 1,
      YCoordinate = issue.range.XY$y_coords[i],#terrSize.max, # in full cpp
      terrSizeMax = terrSize.max
    )
  )
  print(i)
  print(outputs.spread)
}
