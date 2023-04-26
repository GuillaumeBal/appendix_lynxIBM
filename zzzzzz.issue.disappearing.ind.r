# indivivuals disappeared  already from from nextCellsType_whoSorted
# found from  "Yes and No corr move" sanity loop
issue.outputs <- 
  outputs.loop[[matching.size.issue[1]]]


issue.outputs$nDispLeft
issue.outputs$Dispersers_who %>% unique %>% length

issue.outputs$nextCellsType_whoSorted %>% unique %>% length
issue.outputs$nextCellsType_who %>% unique %>% length

issue.outputs$nextCellsType_whoSorted[issue.outputs$nextCellsType_IsMoveCorrF == 1] %>% unique
issue.outputs$nextCellsType_whoSorted[issue.outputs$nextCellsType_IsMoveCorrF == 0] %>% unique
issue.outputs$CellsDisp_who %>% unique %>% length
