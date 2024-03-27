source("./utils/global.R")
source("./utils/enrichment.R")

################ ORA - GO-terms ################################################

### ORA - GO-terms

#' @description
#' Perform whole enrichment analysis on GO terms with ORA
#' @param reactive_annotated_data : data.frame - data.frame from Rshiny reactive filter
#' @param organism_db : organism annotation db to yse (example : "org.Hs.eg.db")
#' @param ontology : ontologies to consider for GO-term analysis. Should use "MF", "CC" and "BP" separately rather than using "ALL"
#' @param p_value_cutoff : p value cutoff
#' @param p_adj_cutoff :  adjusted p value cutoff
#' @param q_value_cutoff :  adjusted q value cutoff
#' @return enrichResults - Ready to plot enrichment analysis results
#'
#' @example do_ora_go_terms(filtered_data,"org.Hs.eg.db", "MF", 0.01, 0.05, 0.05)
#'
do_ora_go_terms <- function(reactive_annotated_data, organism_db, ontology, p_value_cutoff, p_adj_cutoff, q_value_cutoff) {
  
  ora_ids <- prepare_ora(table_filtered_new_ids)
  universe <- prepare_universe(reactive_annotated_data, organism_db, from = "ENSEMBL")
  
  ora_go <- load_ora_go(ora_ids, universe, organism_db)
  ora_go_after_filter <- filter_table_enrich_results(ora_go, 
                                                     p_value_cutoff = p_value_cutoff, 
                                                     p_adj_cutoff   = p_adj_cutoff, 
                                                     q_value_cutoff = q_value_cutoff)
  results <- filter_go_enrich_results(ora_go_after_filter, ontology = ontology)
  
  return(results)
}

#######################

