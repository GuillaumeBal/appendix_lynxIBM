cell.to.check <- sample.int(n = length(cell.index), size = 1000, replace = FALSE)
gb.cells.coord <- CellNumtoRowCol(cellNum = cell.to.check, Matrix = avail.mat.check,
                                  cpp_as_input = TRUE, cpp_as_output = 0) %>% as.data.frame()
gb.cells.coord %>% head
cell.to.check - gb.cells.coord$cellNum

terra.cells.coord <- terra::xyFromCell(cell = cell.to.check + 1, object = sim$availCellsRas) %>% 
  as.data.frame()

checking.xy.df <- data.frame(
  gb = paste(gb.cells.coord$x_coords, gb.cells.coord$y_coords, sep = '_'), 
  terra = paste(terra.cells.coord$x, terra.cells.coord$y, sep = '_')
)

checking.xy.df %>% head(5)

table(terra.cells.coord$y == gb.cells.coord$y_coords)
table(terra.cells.coord$x == gb.cells.coord$x_coords)


checking.xy.df$gb %>% 
 `==` (checking.xy.df$terra) %>% table
