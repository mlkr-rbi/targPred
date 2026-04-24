library(rols)
library(dplyr)
library(stringr)



# Downloading efo ontology 
# 21.02. This is obsolete - using previously generated results
efo_terms <- rols::terms("EFO", pagesize=1000)

#


# Creating pairwise mappings

## Initializing empty df
efo_mappings <- data.frame(id1=NA, id2=NA, label=NA) 

#cbind(efo_mappings, data.frame(id, unlist(efo_terms[[id]]@annotation$database_cross_reference), efo_terms[[id]]@label))

## Looping through elements of the list
for(i in 1:length(names(efo_terms@x))){
  id <- names(efo_terms@x)[i]
  if(i%%100 == 0) {message(i)}
  
  if(is.null(unlist(efo_terms[[id]]@annotation$database_cross_reference))){
    ## Some ids have no cross refrences
    working_df <-  data.frame(id, NA, efo_terms[[id]]@label)
  }
  
  else{
    working_df <-  data.frame(id, unlist(efo_terms[[id]]@annotation$database_cross_reference), efo_terms[[id]]@label)
  }
  
  ## Changing colnames of working_df before binding
  colnames(working_df) <- c("id1","id2","label")
  efo_mappings <- rbind(efo_mappings, working_df)
}

# Remove the initial NA row
efo_mappings <- efo_mappings %>% filter(!is.na(id1))

saveRDS(efo_terms, "data/processed/efo_terms.rds")
write.csv(efo_mappings, "data/processed/efo_mappings.csv", row.names = FALSE)
