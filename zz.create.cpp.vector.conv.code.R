colnames(lynx.gb)

lynx.gb <- sim$lynx[,]
disperser.gb <- lynx.gb[lynx.gb$status == 'disp', ]

# search territory from lynx dataframe to vectors
num.col.lynx <- lynx.gb %>%  `[`(1, ) %>% sapply(., is.numeric) 

# get lynx data into vectors
paste('lynx_', colnames(lynx.gb)[num.col.lynx], ' = lynx["', colnames(lynx.gb)[num.col.lynx], '"]', sep = "") %>% 
  paste(., collapse = ', ') %>% paste('NumericVector ', ., ';', sep = '') %>% 
  paste(., '\n',
        paste('lynx_', colnames(lynx.gb)[!num.col.lynx], ' = lynx["', colnames(lynx.gb)[!num.col.lynx], '"]', sep = "") %>% 
          paste(., collapse = ', ') %>% paste('CharacterVector ', ., ';', sep = '')) %>% cat() 

# define residents and disperser vectors
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

# update through loop ====================================================================
paste('lynx_', colnames(lynx.gb), '(i) = dispersers_', colnames(lynx.gb), '_new(i-residents_xcor.size())', sep = "") %>% 
  paste(., collapse = ', ') %>% paste( ., ';', sep = '') %>% cat

paste('lynx_', colnames(lynx.gb), '(i) = residents_', colnames(lynx.gb),  '(i - dispersers_who_new.size())', sep = "") %>% 
  paste(., collapse = ', ') %>% paste( ., ';', sep = '') %>% cat

# remove lynx "dead" from end of table, in fact reduce the table that's been reorganized ============
paste('lynx_', colnames(lynx.gb), '.erase(nLynx_new, nLynx - 1)', sep = '') %>% 
  paste(., collapse = ', ') %>% paste( ., ';', sep = '') %>% cat

# use of clone for lynx table =======================================================================

# define residents and disperser vectors
paste('lynx_', colnames(lynx.gb)[num.col.lynx], '_new(nLynx_new)', sep = "") %>% paste(., collapse = ', ') %>% # resideents (nLynx - nDisp)
  paste('IntegerVector ', ., ';', sep = '') %>% 
  paste(., '\n',
        paste('lynx_', colnames(lynx.gb)[!num.col.lynx], '_new(nLynx_new)', sep = "") %>% 
          paste(., collapse = ', ') %>% paste('CharacterVector ', ., ';', sep = '')) %>% cat() 


paste('lynx_', colnames(lynx.gb), '_new[res_index] = residents_',  colnames(lynx.gb), sep = "") %>% 
  paste(., collapse = ', ') %>% paste( ., ';', sep = '') %>% cat
lynx_xcor_new[res_index] = residents_xcor;

paste('lynx_', colnames(lynx.gb), '_new[res_index] = residents_',  colnames(lynx.gb), sep = "") %>% 
  paste(., collapse = ', ') %>% paste( ., ';', sep = '') %>% cat

paste('lynx_', ChosenCells_variables, '_new[disp_new_index] = dispersers_',  ChosenCells_variables, "_new", sep = "") %>% 
  paste(., collapse = ', ') %>% paste( ., ';', sep = '') %>% cat

varNotChosenCells <- 
colnames(lynx.gb)[colnames(lynx.gb) %>% `%in%`(ChosenCells_variables) %>% `!`]
paste('lynx_', varNotChosenCells, '_new(i + nRes) = lynx_',  varNotChosenCells, "(corresLynx)", sep = "") %>% 
  paste(., collapse = ', ') %>% paste( ., ';', sep = '') %>% cat


# for dispersers with chosen cells var
lynx_xcor_new[disp_new_index] = dispersers_xcor_new;
# for other var notin chosen 


paste('lynx_', colnames(lynx.gb), '= clone(lynx_', colnames(lynx.gb), '_new)', sep = '') %>% 
  paste(., collapse = ', ') %>% paste( ., ';', sep = '') %>% cat


# export all var to table ======================================================
paste('_["dispersers_', colnames(lynx.gb), '_new"] = dispersers_', colnames(lynx.gb), '_new', sep = "") %>% 
  paste(., collapse = ', ') %>% cat
paste('_["residents_', colnames(lynx.gb), '"] = residents_', colnames(lynx.gb), sep = "") %>% 
  paste(., collapse = ', ') %>% cat

# reshape all dispersre data to smaller size if needed=========================

dispersers_name.erase(nDispOld, nDisp - 1);

paste('dispersers_', c(colnames(lynx.gb), 'steps'), '.erase(nDispOld, nDisp - 1)', sep = "") %>% 
  paste(., collapse = ', ') %>% paste(., ';', sep = '') %>% cat

paste('residents_', c(colnames(lynx.gb), 'steps'), '.resize(nRes)', sep = "") %>% 
  paste(., collapse = ', ') %>% paste(., ';', sep = '') %>% cat


paste('residents_', c(colnames(lynx.gb), 'steps'), '.push_back(lynx_', c(colnames(lynx.gb), 'steps'), '(i))', sep = '') %>% 
  paste(., collapse = ', ') %>% paste( ., ';', sep = '') %>% cat


