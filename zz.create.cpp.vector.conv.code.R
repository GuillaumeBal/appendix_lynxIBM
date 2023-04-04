colnames(disp.gb)

# search territory from Disp dataframe to vectors
num.col.disp <- disp.gb %>%  `[`(1, ) %>% sapply(., is.numeric) 
paste('Disp_', colnames(disp.gb)[num.col.disp], ' = Disp["', colnames(disp.gb)[num.col.disp], '"]', sep = "") %>% 
  paste(., collapse = ', ') %>% paste('NumericVector ', ., ';', sep = '') %>% 
  paste(., '\n',
        paste('Disp_', colnames(disp.gb)[!num.col.disp], ' = Disp["', colnames(disp.gb)[!num.col.disp], '"]', sep = "") %>% 
          paste(., collapse = ', ') %>% paste('CharacterVector ', ., ';', sep = '')) %>% cat() 

# to create Disp_fem subset vectors
num.col.disp <- disp.gb %>%  `[`(1, ) %>% sapply(., is.numeric) 
paste('DispFem_', colnames(disp.gb)[num.col.disp], ' = Disp_', colnames(disp.gb)[num.col.disp], '[index_fem_Disp]', sep = "") %>% 
  paste(., collapse = ', ') %>% paste('NumericVector ', ., ';', sep = '') %>% 
  paste(., '\n',
        paste('DispFem_', colnames(disp.gb)[!num.col.disp], ' = Disp_', colnames(disp.gb)[!num.col.disp], '[index_fem_Disp]', sep = "") %>% 
          paste(., collapse = ', ') %>% paste('CharacterVector ', ., ';', sep = '')) %>% cat() 
