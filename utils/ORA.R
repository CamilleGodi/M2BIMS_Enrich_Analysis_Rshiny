################ ORA - GO-terms ################################################

### ORA - GO-terms

#' @description
#' Perform whole enrichment analysis on GO terms with ORA
#' @param reactive_annotated_data : data.frame - data.frame from Rshiny reactive filter
#' @param organism_db : organism annotation db to yse (example : "org.Hs.eg.db")
#' @param universe : organism "universe" prepared with prepare_universe(reactive_annotated_data, organism_db, "ENSEMBL") 
#' @return enrichResults - Ready to plot enrichment analysis results
#'
#' @example do_ora_go_terms(filtered_data,"org.Hs.eg.db")
#'
do_ora_go_terms <- function(reactive_annotated_data, organism_db, universe) {
  ora_ids <- prepare_ora(reactive_annotated_data)
  ora_go <- load_ora_go(ora_ids, universe, organism_db)
  return(ora_go)
}

#######################

### ORA - KEGG

#' @description
#' Perform whole enrichment analysis on GO terms with ORA
#' @param reactive_annotated_data : data.frame - data.frame from Rshiny reactive filter
#' @param organism_db : organism annotation db to yse (example : "org.Hs.eg.db")
#' @param universe : organism "universe" prepared with prepare_universe(reactive_annotated_data, organism_db, "ENSEMBL") 
#' @param kegg_organism_code : organism 3- nor 4-letters code (example : "hsa")sis results
#'
#' @example do_ora_kegg(filtered_data, "org.Hs.eg.db", "hsa")
#'
do_ora_kegg <- function(reactive_annotated_data, organism_db, universe, kegg_organism_code) {
  ora_ids <- prepare_ora(reactive_annotated_data)
  ora_kegg <- load_ora_kegg(gene_list = ora_ids,
                            organism_db = kegg_organism_code,
                            universe = universe)
  return(ora_kegg)
}

#######################

### ORA - Reactome

#' @description
#' Perform whole enrichment analysis on GO terms with ORA
#' @param reactive_annotated_data : data.frame - data.frame from Rshiny reactive filter
#' @param reactome_organism_name : organism name in reactome (example : "human")
#' @param universe : organism "universe" prepared with prepare_universe(reactive_annotated_data, organism_db, "ENSEMBL") 
#' @return enrichResults - Ready to plot enrichment analysis results
#'
#' @example do_ora_reactome(filtered_data, "human")
#'
do_ora_reactome <- function(reactive_annotated_data, reactome_organism_name, universe) {
  ora_ids <- prepare_ora(reactive_annotated_data)
  ora_reactome <- load_ora_reactome(gene_list = ora_ids,
                                    organism_db = reactome_organism_name,
                                    universe = universe)
  return(ora_reactome)
}