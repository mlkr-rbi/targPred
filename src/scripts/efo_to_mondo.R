
library(DODO)
connect_to_dodo("https://genodesy.org/DODO/", remember=TRUE)

project_dir <- "/home/katarina/MyProjects/targPred"
inputs_dir <- "data/inputs"
processed_dir <- "data/processed"
outputs_dir <- "/data/outputs"

gwas_grb_clnv <- read.csv(file = file.path(project_dir, outputs_dir, "loci_in_GRBs_clnv.csv"))


efo_ids_to_map <- gwas_grb_clnv %>%
  filter(is.na(MONDO_ID)) %>%
  pull(MAPPED_TRAIT_ID) %>%
  unique()


dodo_efo_mondo <- list()
ids_err <- c()

for(i in 1:length(efo_ids_to_map)){
  
  # This looks for direct references
  dodo_efo_mondo[[efo_ids_to_map[i]]] <- convert_concept(from = efo_ids_to_map[i],
                                                         to = "MONDO")
  
  # 
  if(dim(dodo_efo_mondo[[efo_ids_to_map[i]]])[1] == 0){
    dodo_efo_mondo[[efo_ids_to_map[i]]] <- convert_concept(from = efo_ids_to_map[i], 
                                                           to = "MONDO",
                                                           step = NULL,
                                                           intransitive_ambiguity = NULL)
    
  }
  if(dim(dodo_efo_mondo[[efo_ids_to_map[i]]])[1] == 0){
    tryCatch({
      dodo_efo_mondo[[efo_ids_to_map[i]]] <- get_related(from = efo_ids_to_map[i], 
                                                         to = "MONDO",
                                                         from.concept = "Phenotype",
                                                         to.concept = "Phenotype",
                                                         step = NULL,
                                                         intransitive_ambiguity = NULL)
    }, error = function(e) {
      # Handle the error here (you can print a message, skip, or take other action)
      message(paste("Error for", efo_ids_to_map[i], ":", conditionMessage(e)))
      ids_err <- c(ids_err, efo_ids_to_map[i])
    })
  }
  
  if(i%%10 == 0 || i == length(efo_ids_to_map)){
    message(paste("[Processing", i,"out of", length(efo_ids_to_map), "]"))
    saveRDS(dodo_efo_mondo, file = file.path(project_dir, processed_dir, "efo_to_mondo_mappigns_with_DODO.rds"))
    saveRDS(ids_err, file = file.path(project_dir, processed_dir, "efo_to_mondo_failed.rds"))
  }
  
}

