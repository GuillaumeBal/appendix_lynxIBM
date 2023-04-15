colnames(lynx.gb)

# search territory from lynx dataframe to vectors
num.col.lynx <- lynx.gb %>%  `[`(1, ) %>% sapply(., is.numeric) 

# get lynx data into vectors
paste('lynx_', colnames(lynx.gb)[num.col.lynx], ' = lynx["', colnames(lynx.gb)[num.col.lynx], '"]', sep = "") %>% 
  paste(., collapse = ', ') %>% paste('NumericVector ', ., ';', sep = '') %>% 
  paste(., '\n',
        paste('lynx_', colnames(lynx.gb)[!num.col.lynx], ' = lynx["', colnames(lynx.gb)[!num.col.lynx], '"]', sep = "") %>% 
          paste(., collapse = ', ') %>% paste('CharacterVector ', ., ';', sep = '')) %>% cat() 

# define reidents and disperser vectors
paste('dispersers_', colnames(lynx.gb)[num.col.lynx], '(nDisp)', sep = "") %>% paste(., collapse = ', ') %>% # resideents (nLynx - nDisp)
  paste('IntegerVector ', ., ';', sep = '') %>% 
  paste(., '\n',
        paste('dispersers_', colnames(lynx.gb)[!num.col.lynx], '(nDisp)', sep = "") %>% 
          paste(., collapse = ', ') %>% paste('CharacterVector ', ., ';', sep = '')) %>% cat() 

# fill some residents and disperser vector through loops
paste('residents_', colnames(lynx.gb), '(i_ndisp) = lynx_', colnames(lynx.gb), '(i)', sep = "") %>% 
  paste(., collapse = ', ') %>% paste(., ";", sep = '') %>% cat


# search territory from lynx dataframe to vectors
paste('lynx_', colnames(lynx.gb)[num.col.lynx], ' = lynx["', colnames(lynx.gb)[num.col.lynx], '"]', sep = "") %>% 
  paste(., collapse = ', ') %>% paste('NumericVector ', ., ';', sep = '') %>% 
  paste(., '\n',
        paste('lynx_', colnames(lynx.gb)[!num.col.lynx], ' = lynx["', colnames(lynx.gb)[!num.col.lynx], '"]', sep = "") %>% 
          paste(., collapse = ', ') %>% paste('CharacterVector ', ., ';', sep = '')) %>% cat() 


# to create disperser_fem subset vectors
num.col.disperser <- disperser.gb %>%  `[`(1, ) %>% sapply(., is.numeric) 
paste('disperserFem_', colnames(disperser.gb)[num.col.disperser], ' = disperser_', colnames(disperser.gb)[num.col.disperser], '[index_fem_disperser_rand]', sep = "") %>% 
  paste(., collapse = ', ') %>% paste('NumericVector ', ., ';', sep = '') %>% 
  paste(., '\n',
        paste('disperserFem_', colnames(disperser.gb)[!num.col.disperser], ' = disperser_', colnames(disperser.gb)[!num.col.disperser], '[index_fem_disperser_rand]', sep = "") %>% 
          paste(., collapse = ', ') %>% paste('CharacterVector ', ., ';', sep = '')) %>% cat() 

# create dataframe from vectors
paste(
  paste(
    'DataFrame disperserFem = DataFrame::create(Named("', colnames(disperser.gb)[1], '") = disperserFem_', colnames(disperser.gb)[1], ', ', sep = ''),
  paste(
    '_["',colnames(disperser.gb)[2:length(colnames(disperser.gb))], '"] = disperserFem_', colnames(disperser.gb)[2:length(colnames(disperser.gb))], sep = '') %>% 
    paste(., collapse = ', '),
  ');', sep = '') %>% cat


# create dispersers new table ===========================================================

ChosenCells_variables <- c(
  'who',
  'ind',
  'hab',
  'pxcorHere',
  'pycorHere',
  'pxcor',
  'pycor',
  'lastDispX',
  'lastDispY',
  'nMat')

colnames.not.in.chosenCells <-
  colnames(lynx.gb) %>% '%in%'(ChosenCells_variables) %>% `!`

paste('dispersers_', ChosenCells_variables, '_new = ChosenCells_',
      ChosenCells_variables, '[pos_alive_who_ord]', sep = '') %>% 
  paste(., collapse = ', ') %>% paste0(., ';') %>% cat

paste('dispersers_', colnames(lynx.gb)[num.col.lynx & colnames.not.in.chosenCells], '_new = dispersers_', colnames(lynx.gb)[num.col.lynx & colnames.not.in.chosenCells], '[index_dispersers_dispersers_new]', sep = "") %>% 
  paste(., collapse = ', ') %>% paste('NumericVector ', ., ';', sep = '') %>% 
  paste(., '\n',
        paste('dispersers_', colnames(lynx.gb)[!num.col.lynx  & colnames.not.in.chosenCells], '_new = dispersers_', colnames(lynx.gb)[!num.col.lynx & colnames.not.in.chosenCells], '[index_dispersers_dispersers_new]', sep = "") %>% 
          paste(., collapse = ', ') %>% paste('CharacterVector ', ., ';', sep = '')) %>% cat() 

# update lynx variables =========================================================
paste('lynx_', colnames(lynx.gb), '[disp_new_index] = dispersers_', colnames(lynx.gb), '_new', sep = "") %>% 
  paste(., collapse = ', ') %>% paste( ., ';', sep = '') %>% cat

paste('lynx_', colnames(lynx.gb), '[res_index] = residents_', colnames(lynx.gb), sep = "") %>% 
  paste(., collapse = ', ') %>% paste( ., ';', sep = '') %>% cat

# remove lynx "dead" from end of table, in fact reduce the table that's been reorganized
paste('lynx_', colnames(lynx.gb), '.erase(nLynx_new, nLynxv- 1)') %>% 
  paste(., collapse = ', ') %>% paste( ., ';', sep = '') %>% cat


# export all var to table ======================================================
paste('_["dispersers_', colnames(lynx.gb), '_new"] = dispersers_', colnames(lynx.gb), '_new', sep = "") %>% 
  paste(., collapse = ', ') %>% cat
paste('_["residents_', colnames(lynx.gb), '"] = residents_', colnames(lynx.gb), sep = "") %>% 
  paste(., collapse = ', ') %>% cat








